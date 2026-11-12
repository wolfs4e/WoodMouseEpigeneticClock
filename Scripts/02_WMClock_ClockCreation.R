################################################################################
# 02_WMClock_ClockCreation.R
# Wood mouse epigenetic clock | elastic-net training and external validation
################################################################################
# PURPOSE
#   Train an age-prediction clock from laboratory wood mouse ear DNA methylation,
#   evaluate it on a predefined test set, compare predictions with three
#   pan-mammalian clocks, and assess transferability across tissue types.
#
# RUN FROM THE R PROJECT ROOT (WoodMouseEpigeneticClock/)
# INPUTS
#   Data/Processed Data/beta_pheno_anno.Rdata  (object: joined2mcpg)
#   Data/PanMamManifest_Aligned_ASylvaticus.csv
#   Data/weights_panmamC1.csv
#   Data/weights_panmamC2.csv
#   Data/weights_panmamC3.csv
#
# OUTPUTS
#   Data/Processed Data/WMClock_cpg_coefficients.Rdata
#   Data/Processed Data/WMClock_cpg_coefficients.csv
#   Data/Processed Data/WMClock_fitted_model.rds (optional reusable glmnet fit)
#   Data/Processed Data/WMClock_Figure1_plots.rds (Figure 1A-C plot objects)
#   Results/02_WMClock_ClockCreation/Tables/  (sample summary, 4-sheet XLSX)
#   Results/02_WMClock_ClockCreation/Figures/ (original-style PDF figures)
#
# METHOD NOTES
#   * Training/test membership is supplied by the metadata column 'Set'.
#   * CpG features are columns beginning 'cg'. Alpha is selected automatically.
#   * Alpha and lambda are selected using leave-one-sample-out cross-validation
#     by minimising cross-validated mean squared prediction error.
#   * Training-fit error is IN-SAMPLE; test accuracy is the validation metric.
#   * Figure 1A-C, Figure S1 and Figure S2 are generated here; S3 is not.
#   * Pan-mammalian clocks use available CpGs and their published transforms.
#   * Tissue comparisons are exploratory applications of the ear-trained clock.
################################################################################

# ---- 1. Packages and project paths -------------------------------------------
suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(openxlsx)
})

# Run with the RStudio project / repository root as the working directory.
input_data <- file.path("Data", "Processed Data", "beta_pheno_anno.Rdata")
manifest_file <- file.path("Data", "PanMamManifest_Aligned_ASylvaticus.csv")
weights_files <- file.path("Data", paste0("weights_panmamC", 1:3, ".csv"))
processed_dir <- file.path("Data", "Processed Data")
result_dir <- file.path("Results", "02_WMClock_ClockCreation")
table_dir <- file.path(result_dir, "Tables")
figure_dir <- file.path(result_dir, "Figures")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(input_data, manifest_file, weights_files)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

# Fix the random cross-validation fold assignment for reproducibility.
set.seed(2026)

# ---- 2. Import methylation data and define sample sets ------------------------
load(input_data)  # loads 'joined2mcpg': samples (rows), metadata and CpGs (columns)
stopifnot(exists("joined2mcpg"), is.data.frame(joined2mcpg))
required_cols <- c("Set", "Age.Days", "Sex", "Mouse.ID", "Batch", "Tissue", "Where")
if (!all(required_cols %in% names(joined2mcpg))) {
  stop("Missing expected columns: ",
       paste(setdiff(required_cols, names(joined2mcpg)), collapse = ", "))
}

# Summarise the predefined training and test groups.
set_summary <- joined2mcpg %>%
  group_by(Set) %>%
  summarise(n = n(),
            mean_age = mean(Age.Days, na.rm = TRUE),
            se_age = sd(Age.Days, na.rm = TRUE) / sqrt(sum(!is.na(Age.Days))),
            min_age = min(Age.Days, na.rm = TRUE),
            max_age = max(Age.Days, na.rm = TRUE),
            n_males = sum(Sex == "m", na.rm = TRUE),
            n_females = sum(Sex == "f", na.rm = TRUE),
            .groups = "drop")
write.csv(set_summary, file.path(table_dir, "WMClock_sample_set_summary.csv"), row.names = FALSE)

train <- joined2mcpg %>% filter(Set == "train")
test <- joined2mcpg %>% filter(Set == "test")
if (!nrow(train) || !nrow(test)) stop("Training or test set is empty; check 'Set'.")

cpg_cols <- grep("^cg", names(train), value = TRUE)
if (!length(cpg_cols)) stop("No CpG columns beginning with 'cg' found.")
if (anyNA(train[, cpg_cols]) || anyNA(train$Age.Days)) {
  stop("Training CpGs and ages must be nonmissing for glmnet.")
}
x_mat_all <- as.matrix(train[, cpg_cols, drop = FALSE])
storage.mode(x_mat_all) <- "double"
y_vec <- train$Age.Days

# ---- 3. Fit elastic-net models across mixing parameters ----------------------
# alpha = 0.1 approaches ridge regression; alpha = 1 is lasso regression.
# Each alpha uses leave-one-sample-out CV to select lambda.min.
# Compare alphas using the CV mean squared error at their selected lambda;
# unlike in-sample MdAE, this evaluates predictions for held-out samples.
alpha_grid <- seq(0.1, 1, by = 0.1)
fit_results <- vector("list", length(alpha_grid))
for (i in seq_along(alpha_grid)) {
  a <- alpha_grid[i]
  message("Fitting wood mouse clock at alpha = ", a)
  cv_fit <- cv.glmnet(x = x_mat_all, y = y_vec, alpha = a,
                      nfolds = nrow(train), grouped = FALSE)
  gl_fit <- glmnet(x = x_mat_all, y = y_vec, alpha = a,
                   lambda = cv_fit$lambda.min)
  raw_coefs <- coef(gl_fit)
  coef_df <- data.frame(CpG = rownames(raw_coefs),
                        Coefficient = as.numeric(raw_coefs)) %>%
    filter(Coefficient != 0)
  intercept_i <- coef_df$Coefficient[coef_df$CpG == "(Intercept)"]
  cg_coef <- coef_df %>% filter(CpG != "(Intercept)")
  if (!nrow(cg_coef)) {
    pred_train <- rep(intercept_i, length(y_vec))
  } else {
    pred_train <- as.numeric(x_mat_all[, cg_coef$CpG, drop = FALSE] %*%
                               cg_coef$Coefficient + intercept_i)
  }
  cv_lambda_index <- which.min(abs(cv_fit$lambda - cv_fit$lambda.min))
  cv_mse <- cv_fit$cvm[cv_lambda_index]
  fit_results[[i]] <- list(alpha = a, n_cpg = nrow(cg_coef),
                           cv_mse = cv_mse,
                           r = if (sd(pred_train) > 0) cor(pred_train, y_vec) else NA_real_,
                           mdae = median(abs(pred_train - y_vec)),
                           coef_table = coef_df, model = gl_fit,
                           cv_model = cv_fit)
}

alpha_summary <- bind_rows(lapply(fit_results, function(z) {
  data.frame(alpha = z$alpha, n_CpGs = z$n_cpg,
             LOOCV_MSE = z$cv_mse,
             training_Pearson_r = z$r, training_MdAE_days = z$mdae)
}))
print(alpha_summary)

# Diagnostic plots use the original line/point styling and minimalist theme.
alpha_error_plot <- ggplot(alpha_summary, aes(x = alpha, y = training_MdAE_days)) +
  geom_line() + geom_point() +
  labs(title = "Median Absolute Error vs Alpha", x = "Alpha (Elastic Net Mixing)",
       y = "Median Absolute Error") + theme_minimal()
alpha_r_plot <- ggplot(alpha_summary, aes(x = alpha, y = training_Pearson_r)) +
  geom_line() + geom_point() + theme_minimal() +
  labs(title = "Correlation vs Alpha", x = "Alpha",
       y = "Pearson Correlation (Epigenetic Age ~ Chronological Age)")
pdf(file.path(figure_dir, "WMClock_alpha_diagnostics.pdf"), width = 9, height = 4)
print(ggarrange(alpha_error_plot, alpha_r_plot, ncol = 2))
dev.off()

# Select the alpha with the smallest held-out CV mean squared prediction error.
# In a tie, which.min() selects the first (smallest) alpha in alpha_grid.
if (any(!is.finite(alpha_summary$LOOCV_MSE))) {
  stop("Non-finite LOOCV error encountered; cannot select alpha.")
}
best_alpha_index <- which.min(alpha_summary$LOOCV_MSE)
selected_fit <- fit_results[[best_alpha_index]]
message(sprintf("Selected alpha = %.1f (LOOCV MSE = %.4f; lambda = %.6g; %d CpGs)",
                selected_fit$alpha, selected_fit$cv_mse,
                selected_fit$cv_model$lambda.min, selected_fit$n_cpg))
best_coefs <- selected_fit$coef_table
saveRDS(selected_fit$model, file.path(processed_dir, "WMClock_fitted_model.rds"))

# ---- 4. Annotate and export clock coefficients -------------------------------
# Preserve the intercept as a separate coefficient row.
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"CGid" %in% names(manifest)) stop("Manifest must contain column 'CGid'.")
WMClock_cpg_coefficients <- left_join(best_coefs, manifest, by = c("CpG" = "CGid"))
if ("annotation2" %in% names(WMClock_cpg_coefficients)) {
  WMClock_cpg_coefficients$annotation2 <- factor(
    WMClock_cpg_coefficients$annotation2,
    levels = c("Downstream", "3' UTR", "Intron", "Exon", "5' UTR",
               "Promoter", "Intergenic"))
}
save(WMClock_cpg_coefficients,
     file = file.path(processed_dir, "WMClock_cpg_coefficients.Rdata"))
write.csv(WMClock_cpg_coefficients,
          file.path(processed_dir, "WMClock_cpg_coefficients.csv"), row.names = FALSE)

# Apply the clock's coefficient vector to any dataset with the needed CpGs.
# Stop rather than silently calculate incomplete predictions if any sites are absent.
predict_wmclock <- function(data, coefs) {
  sites <- setdiff(coefs$CpG, "(Intercept)")
  missing <- setdiff(sites, names(data))
  if (length(missing)) stop("Missing wood mouse clock CpGs: ", paste(missing, collapse = ", "))
  intercept <- coefs$Coefficient[coefs$CpG == "(Intercept)"]
  if (length(intercept) != 1L) stop("Expected exactly one clock intercept.")
  if (!length(sites)) return(rep(intercept, nrow(data)))
  as.numeric(as.matrix(data[, sites, drop = FALSE]) %*%
               coefs$Coefficient[match(sites, coefs$CpG)] + intercept)
}

# Small helper for accuracy summaries. Use the observed, not hard-coded, metrics.
performance <- function(data, predicted_col, label) {
  ok <- is.finite(data[[predicted_col]]) & is.finite(data$Age.Days)
  x <- data[[predicted_col]][ok]
  y <- data$Age.Days[ok]
  data.frame(Dataset = label, n = length(x),
             MdAE = if (length(x)) median(abs(x - y)) else NA_real_,
             Pearson_r = if (length(x) >= 3 && sd(x) > 0 && sd(y) > 0)
               cor(x, y) else NA_real_)
}

# Reproduce the manuscript-style panels: white background with gridlines,
# open jittered circles, solid black regression line and dashed identity line.
# Only the annotation values are computed automatically from each plot's data.
clock_plot <- function(data, predicted_col, label = NULL,
                       x_limits = c(50, 400), y_limits = c(50, 400),
                       y_breaks = seq(100, 400, 100), annotation_y = 400) {
  stats <- performance(data, predicted_col, if (is.null(label)) "" else label)
  annotation <- sprintf("MdAE: %.2f days | r: %.2f", stats$MdAE, stats$Pearson_r)
  ggplot(data, aes(x = Age.Days, y = .data[[predicted_col]])) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_jitter(pch = 1, size = 4,
                position = position_jitter(width = 0.1, height = 0.1)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(x = "Chronological Age (Days)", y = "Epigenetic Age (Days)") +
    scale_x_continuous(breaks = seq(100, 400, 100), limits = x_limits) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    annotate("text", x = 225, y = annotation_y,
             label = annotation, color = "black", size = 3)
}

# ---- 5. Evaluate wood mouse clock on training and test samples ---------------
train$predicted_age_best_model <- predict_wmclock(train, best_coefs)
test$predAllAge <- predict_wmclock(test, best_coefs)

# The metadata's 'A'/'B' identifiers flag repeated samples in this dataset.
# This subset follows that naming convention; verify individual-level
# independence against the underlying mouse ID records before publication.
testsub <- test %>% filter(!grepl("A|B", Mouse.ID))
if (!nrow(testsub)) warning("No independent-identifier test samples remain.")

wm_performance <- bind_rows(
  performance(train, "predicted_age_best_model", "Training (in-sample)"),
  performance(test, "predAllAge", "Full held-out test set"),
  performance(testsub, "predAllAge", "Test subset excluding A/B identifiers")
)
print(wm_performance)

# Confirm whether the subset contains any IDs also found in the training data.
overlap_ids <- intersect(unique(as.character(testsub$Mouse.ID)),
                         unique(as.character(train$Mouse.ID)))
if (length(overlap_ids)) {
  warning("Test subset includes Mouse.ID(s) present in training: ",
          paste(overlap_ids, collapse = ", "))
}

All_graph <- clock_plot(train, "predicted_age_best_model", "Training")
TestAll <- clock_plot(test, "predAllAge", "Test")
if (nrow(testsub) >= 3) {
  TestSub <- clock_plot(testsub, "predAllAge", "Independent test")
}
# Figure 1 includes additional panels (D-G) produced by other analyses.
# This script generates the clock-performance panels (A-C) only.

# ---- 6. Assess availability of pan-mammalian clock CpGs ----------------------
# A single workbook contains an overall retention summary and missing CpG IDs
# for each of the three pan-mammalian clocks.
weight_tables <- lapply(weights_files, read.csv, stringsAsFactors = FALSE)
weight_cpg_cols <- c("coef.name", "var", "var")
for (i in seq_along(weight_tables)) {
  if (!weight_cpg_cols[i] %in% names(weight_tables[[i]])) {
    stop("Missing CpG ID column in ", weights_files[i])
  }
}
get_clock_sites <- function(weights, id_col) {
  ids <- unique(as.character(weights[[id_col]]))
  ids[!is.na(ids) & grepl("^cg", ids)]
}
clock_sites <- Map(get_clock_sites, weight_tables, weight_cpg_cols)
missing_sites <- lapply(clock_sites, setdiff, y = names(joined2mcpg))
retention <- bind_rows(lapply(seq_along(clock_sites), function(i) {
  needed <- length(clock_sites[[i]])
  retained <- needed - length(missing_sites[[i]])
  data.frame(Clock = paste0("Pan-mammalian clock ", i),
             n_clock_CpGs = needed, n_retained = retained,
             n_missing = length(missing_sites[[i]]),
             percent_retained = if (needed) 100 * retained / needed else NA_real_)
}))
print(retention)

wb <- createWorkbook()
addWorksheet(wb, "Probe retention")
writeData(wb, "Probe retention", retention)
freezePane(wb, "Probe retention", firstRow = TRUE)
setColWidths(wb, "Probe retention", cols = 1:ncol(retention), widths = "auto")
for (i in seq_along(missing_sites)) {
  sheet <- paste0("Clock ", i, " missing CpGs")
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data.frame(CpG = missing_sites[[i]]))
  freezePane(wb, sheet, firstRow = TRUE)
  setColWidths(wb, sheet, cols = 1, widths = 22)
}
saveWorkbook(wb, file.path(table_dir, "PanMam_clock_CpG_retention.xlsx"), overwrite = TRUE)

# ---- 7. Apply pan-mammalian clocks -------------------------------------------
# Retain the supplied age-transform constants and coefficient ordering.
# Check the retained fraction above before interpreting external clock error.
score_pan <- function(data, weights, id_col, coef_col) {
  if (!coef_col %in% names(weights)) stop("Missing coefficient column: ", coef_col)
  sites <- intersect(as.character(weights[[id_col]]), names(data))
  sites <- unique(sites[grepl("^cg", sites)])
  if (!length(sites)) stop("No overlapping pan-mammalian clock CpGs.")
  beta <- as.matrix(data[, sites, drop = FALSE])
  coefs <- weights[[coef_col]][match(sites, weights[[id_col]])]
  if (anyNA(coefs)) stop("Missing coefficients for retained clock sites.")
  as.numeric(beta %*% coefs + weights[[coef_col]][1])
}

# Published pan-mammalian clock transformations, with values in years.
Lifespan <- 6.3
SexualMaturity <- 0.186301
Gestation <- 0.0630137

pm1_linear <- score_pan(joined2mcpg, weight_tables[[1]], "coef.name", "coef.value")
joined2mcpg$predPM1 <- (exp(pm1_linear) - 2) * 365

pm2_linear <- score_pan(joined2mcpg, weight_tables[[2]], "var", "beta_clock2")
joined2mcpg$predPM2 <-
  (exp(-exp(-pm2_linear)) * (Lifespan + Gestation) - Gestation) * 365

pm3_linear <- score_pan(joined2mcpg, weight_tables[[3]], "var", "beta_clock3")
adultadjusted <- SexualMaturity + Gestation
mhat <- 5 * (Gestation / SexualMaturity)^0.38
joined2mcpg$predPM3 <- ifelse(
  pm3_linear >= 0,
  mhat * adultadjusted * (pm3_linear + 1) - Gestation,
  mhat * adultadjusted * exp(pm3_linear) - Gestation
) * 365

pan_performance <- bind_rows(lapply(1:3, function(i) {
  performance(joined2mcpg, paste0("predPM", i),
              paste0("Pan-mammalian clock ", i, ": all samples"))
}))
print(pan_performance)

lab_colony <- joined2mcpg %>% filter(Where == "Lab", is.finite(Age.Days))
PanOnly_Graph <- clock_plot(lab_colony, "predPM1", "Pan-mammalian clock 1",
                            y_limits = c(0, 1200), y_breaks = seq(100, 1200, 200),
                            annotation_y = 1200)
PanOnly_Graph2 <- clock_plot(lab_colony, "predPM2", "Pan-mammalian clock 2",
                             y_limits = c(0, 1200), y_breaks = seq(100, 1200, 200),
                             annotation_y = 1200)
PanOnly_Graph3 <- clock_plot(lab_colony, "predPM3", "Pan-mammalian clock 3",
                             y_limits = c(0, 1200), y_breaks = seq(100, 1200, 200),
                             annotation_y = 1200)
# Figure 1A-C: this is the three-panel clock comparison, not the full Figure 1.
pdf(file.path(figure_dir, "Figure_1_panels_A-C.pdf"), width = 8, height = 2.66)
print(ggarrange(All_graph, TestAll, PanOnly_Graph2,
                labels = c("A", "B", "C"), ncol = 3))
dev.off()

# Figure S2: three pan-mammalian clocks for known-age colony mice.
pdf(file.path(figure_dir, "Figure_S2.pdf"), width = 8, height = 2.66)
print(ggarrange(PanOnly_Graph, PanOnly_Graph2, PanOnly_Graph3,
                labels = c("A", "B", "C"), ncol = 3))
dev.off()

# Apply pan-mammalian clock 1 to the independent-identifier laboratory subset.
lab_testsub_pm1 <- joined2mcpg %>%
  filter(Where == "Lab", Set == "test", !grepl("A|B", Mouse.ID))
if (nrow(lab_testsub_pm1) >= 3 && nrow(testsub) >= 3) {
  pm1_subset_perf <- performance(lab_testsub_pm1, "predPM1",
                                  "Pan-mammalian clock 1: lab test subset")
  print(pm1_subset_perf)
  PanOnly_GraphTest <- clock_plot(lab_testsub_pm1, "predPM1",
                                   "Pan-mammalian clock 1: lab test subset",
                                   y_limits = c(0, 1200),
                                   y_breaks = seq(100, 1200, 200),
                                   annotation_y = 1200)
  # Figure S1: same two-panel arrangement and plot geometry as manuscript.
  pdf(file.path(figure_dir, "Figure_S1.pdf"), width = 6.91, height = 2.66)
  print(ggarrange(TestSub, PanOnly_GraphTest, labels = c("A", "B")))
  dev.off()
}

# ---- 8. Transfer the ear-trained clock to other tissues -----------------------
# Identify tissue-comparison samples using metadata Batch == 'Horvath'.
tdata <- joined2mcpg %>% filter(Batch == "Horvath")
if (nrow(tdata)) {
  tdata$MouseNum <- as.numeric(gsub("[^0-9]", "", as.character(tdata$Mouse.ID)))
  tdata$predAllAge <- predict_wmclock(tdata, best_coefs)

  tissue_performance <- bind_rows(lapply(split(tdata, tdata$Tissue), function(df) {
    performance(df, "predAllAge", as.character(df$Tissue[1]))
  }))
  names(tissue_performance)[names(tissue_performance) == "Dataset"] <- "Tissue"
  print(tissue_performance)

  TestTissue <- ggplot(tdata, aes(Age.Days, predAllAge)) +
    geom_jitter(aes(color = Tissue), pch = 1, size = 4,
                position = position_jitter(width = 0.1, height = 0.1)) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(x = "Chronological Age (Days)", y = "Epigenetic Age (Days)") +
    scale_x_continuous(breaks = seq(100, 400, 100), limits = c(50, 400)) +
    scale_y_continuous(breaks = seq(100, 400, 100), limits = c(50, 400)) +
    facet_wrap(~Tissue) + theme(legend.position = "none")
  ggsave(file.path(figure_dir, "WMClock_tissue_accuracy.pdf"),
         TestTissue, width = 9, height = 6)

  # Average replicate values per individual/tissue before tissue-pair comparisons.
  df_wide <- tdata %>%
    filter(!is.na(Tissue), !is.na(MouseNum)) %>%
    group_by(MouseNum, Tissue) %>%
    summarise(predWMAge = mean(predAllAge, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Tissue, values_from = predWMAge)
  tissue_names <- setdiff(names(df_wide), "MouseNum")
  if (length(tissue_names) >= 2) {
    tissue_pairs <- combn(tissue_names, 2, simplify = FALSE)
    cor_results <- bind_rows(lapply(tissue_pairs, function(pair) {
      x <- df_wide[[pair[1]]]; y <- df_wide[[pair[2]]]
      valid <- is.finite(x) & is.finite(y)
      n <- sum(valid)
      data.frame(tissue1 = pair[1], tissue2 = pair[2], n = n,
                 Pearson_r = if (n >= 3 && sd(x[valid]) > 0 && sd(y[valid]) > 0)
                   cor(x[valid], y[valid]) else NA_real_)
    }))
    print(cor_results)
  }

  needed_tissues <- c("Ear", "Spleen", "Blood", "Liver")
  if (all(needed_tissues %in% names(df_wide))) {
    plot_data <- df_wide %>%
      select(MouseNum, all_of(needed_tissues)) %>%
      pivot_longer(cols = c(Spleen, Blood, Liver),
                   names_to = "tissue", values_to = "value")
    tissue_compare <- ggplot(plot_data, aes(Ear, value)) +
      geom_point(shape = 1, size = 2.5) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      facet_wrap(~tissue, scales = "free") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      theme_bw() +
      labs(x = "Ear Epigenetic Age (days)",
           y = "Other Tissue Epigenetic Age (days)") +
      scale_x_continuous(breaks = seq(100, 400, 100), limits = c(50, 400)) +
      scale_y_continuous(breaks = seq(100, 400, 100), limits = c(50, 400))
    # Supplemental manuscript assigns S5 to both this and a parasite figure.
    # Use a descriptive name until figure numbering is resolved.
    ggsave(file.path(figure_dir, "WMClock_ear_vs_other_tissues.pdf"),
           tissue_compare, width = 9, height = 4)
  } else {
    warning("Ear/spleen/blood/liver comparison omitted: not all tissue columns exist.")
  }
} else {
  warning("No rows with Batch == 'Horvath'; tissue-transfer analysis skipped.")
}


# ---- 9. Save Figure 1A-C plot objects for assembly with EWAS/GREAT panels -----
# Save native ggplot objects, not screenshots or exported PDF pages. This allows
# Script 03 to assemble the complete manuscript figure directly with ggarrange()
# and export all panels as editable-resolution vector graphics in one PDF.
# The plots retain the data, themes, annotations, and layer settings defined above.
figure1_clock_plots <- list(
  All_graph = All_graph,             # Panel A: training-set clock performance
  TestAll = TestAll,                 # Panel B: held-out test-set performance
  PanOnly_Graph2 = PanOnly_Graph2    # Panel C: pan-mammalian clock 2
)
figure1_plot_file <- file.path(processed_dir, "WMClock_Figure1_plots.rds")
saveRDS(figure1_clock_plots, figure1_plot_file)
message("Saved reusable Figure 1A-C ggplot objects to: ", figure1_plot_file)

message("Clock creation and validation completed.\n",
        "Coefficients: ", processed_dir, "\n",
        "Tables: ", table_dir, "\n",
        "Figures: ", figure_dir)
