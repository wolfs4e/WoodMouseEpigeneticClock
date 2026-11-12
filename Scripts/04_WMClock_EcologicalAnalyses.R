################################################################################
# 04_WMClock_EcologicalAnalyses.R
#
# PURPOSE
#   Apply wood mouse and pan-mammalian epigenetic clocks to field samples and
#   investigate environmental treatments, parasites, phenotypes, population
#   abundance, age structure, and longitudinal epigenetic ageing.
#
# RUN FROM
#   Root of the WoodMouseEpigeneticClock RStudio project.
#
# REQUIRED INPUTS
#   Data/Processed Data/Field20232024_Compiled_wDNAm.Rdata (object: joined)
#   Data/Processed Data/WMClock_cpg_coefficients.csv
#   Data/weights_panmamC1.csv
#   Data/weights_panmamC2.csv
#   Data/weights_panmamC3.csv
#
# OUTPUTS
#   Data/Processed Data/WMClock_EcologicalData.Rdata (object: data)
#   Data/Processed Data/WMClock_LongitudinalData.Rdata
#       (objects: WildLong, WildLong2, WildLong_no_extreme, WildLong_common)
#   Results/04_WMClock_EcologicalAnalyses/Tables/  - sensitivity CSV tables
#   Results/04_WMClock_EcologicalAnalyses/Figures/ - manuscript figures only
#       Figure_2.pdf, Figure_3.pdf, Figure_S3.pdf, Figure_S4.pdf,
#       Figure_S5.pdf (parasite treatment)
#
# ANALYTICAL NOTE
#   Core models, transformations, and original ggplot styling are retained.
#   Original manuscript ggplot styling is retained. The H. polygyrus model
#   uses the convergent zero-inflated Poisson specification and emmeans
#   treatment-by-period contrasts on the expected-response scale.
#   R objects are saved to support later analyses.
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(lme4)
  library(plyr)
  library(dplyr)
  library(purrr)
  library(lmodel2)
  library(optimx)
  library(emmeans)
  library(DHARMa)
  library(pscl)
  library(glmmTMB)
  library(tidyr)
  library(stringr)
  library(visreg)
  library(ggpubr)
  library(Rmisc)
  library(lubridate)
})

# Output folders are created automatically.
processed_dir <- file.path("Data", "Processed Data")
results_dir <- file.path("Results", "04_WMClock_EcologicalAnalyses")
tables_dir <- file.path(results_dir, "Tables")
figures_dir <- file.path(results_dir, "Figures")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
## Formatting the Data #########################################################
################################################################################
# ---- Load full dataset -------------------------------------------------------
load(file.path(processed_dir, "Field20232024_Compiled_wDNAm.Rdata"))
stopifnot(exists("joined"), is.data.frame(joined))
# Ensure subsequent date comparisons, sorting and time intervals use Date objects.
joined$Date <- as.Date(joined$Date)
stopifnot(all(c("EarID", "Date", "AnalysisID", "cg00013189") %in% names(joined)))

# ---- Keep phenotype-only data ------------------------------------------------
pheno <- joined[, !grepl("cg", names(joined))]

# ---- Subset cg columns + EarID, remove rows with missing data in a representative CpG
dnam <- joined[!is.na(joined$cg00013189), grepl("^cg", names(joined)) | names(joined) == "EarID" | names(joined) == "Date"]

# ---- WM Clock ----------------------------------------------------------------
All_weights <- read.csv(file.path(processed_dir, "WMClock_cpg_coefficients.csv"))

a <- intersect(All_weights$CpG, names(dnam)) # order by All_weights$var, not names(WM_data)
data_m_2 <- as.matrix(dnam[, a]) # make matrix of the above intersect
betas <- All_weights$Coefficient[match(a, All_weights$CpG)] #ensuring coefficients align
pred <- data_m_2 %*% betas + All_weights$Coefficient[1] # maxtrix manipulation and adding intercept to each predicted value
dnam$predWMAge <- as.numeric(pred)
hist(dnam$predWMAge)

# ---- PM Clock 1 --------------------------------------------------------------
All_weights <- read.csv(file.path("Data", "weights_panmamC1.csv"))

a <- intersect(All_weights$coef.name, names(dnam)) # order by All_weights$var, not names(WM_data)
data_m_2 <- as.matrix(dnam[, a]) # make matrix of the above intersect
betas <- All_weights$coef.value[match(a, All_weights$coef.name)] #ensuring coefficients align
pred <- data_m_2 %*% betas + All_weights$coef.value[1] # maxtrix manipulation and adding intercept to each predicted value
dnam$predPMAge1 <- as.numeric(pred)
dnam$predPMAge1 <- (exp(dnam$predPMAge1)-2)*365
hist(dnam$predPMAge1)

PMGraph1 <- ggplot(dnam, aes(x = predWMAge, y = predPMAge1)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "WM Clock EpiAge (days)", y = "PM Clock 1 EpiAge (days)")
PMGraph1

# ---- PM Clock 2 --------------------------------------------------------------
Lifespan <-  6.3 # ridiculous for wild mice, closer to 0.5 in wild
SexualMaturity <- 0.186301
Gestation <- 0.0630137

All_weights <- read.csv(file.path("Data", "weights_panmamC2.csv"))

a <- intersect(All_weights$var, names(dnam)) # order by All_weights$var, not names(WM_data)
data_m_2 <- as.matrix(dnam[, a]) # make matrix of the above intersect
betas <- All_weights$beta_clock2[match(a, All_weights$var)] #ensuring coefficients align
pred <- data_m_2 %*% betas + All_weights$beta_clock2[1] # maxtrix manipulation and adding intercept to each predicted value
dnam$predPMAge2 <- as.numeric(pred)

dnam$predPMAge2 <- (exp(-exp(-dnam$predPMAge2))*(Lifespan + Gestation) - Gestation)*365
hist(dnam$predPMAge2)

PMGraph2 <- ggplot(dnam, aes(x = predWMAge, y = predPMAge2)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "WM Clock EpiAge (days)", y = "PM Clock 2 EpiAge (days)")
PMGraph2

# ---- PM Clock 3 --------------------------------------------------------------
Lifespan <-  6.3 # ridiculous for wild mice, closer to 0.5 in wild
SexualMaturity <- 0.186301
Gestation <- 0.0630137

All_weights <- read.csv(file.path("Data", "weights_panmamC3.csv"))

a <- intersect(All_weights$var, names(dnam)) # order by All_weights$var, not names(WM_data)
data_m_2 <- as.matrix(dnam[, a]) # make matrix of the above intersect
betas <- All_weights$beta_clock3[match(a, All_weights$var)] #ensuring coefficients align
pred <- data_m_2 %*% betas + All_weights$beta_clock3[1] # maxtrix manipulation and adding intercept to each predicted value
dnam$predPMAge3 <- as.numeric(pred)

adultadjusted <- SexualMaturity + Gestation
mhat <- 5 * (Gestation/SexualMaturity)^0.38
dnam$predPMAge3 <- (ifelse(dnam$predPMAge3 >= 0, 
                                mhat * adultadjusted * (dnam$predPMAge3 + 1) - Gestation,
                                mhat * adultadjusted * exp(dnam$predPMAge3)  - Gestation))*365
hist(dnam$predPMAge3)

PMGraph3 <- ggplot(dnam, aes(x = predWMAge, y = predPMAge3)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "WM Clock EpiAge (days)", y = "PM Clock 3 EpiAge (days)")
PMGraph3

grDevices::pdf(file.path(figures_dir, "Figure_S3.pdf"), width = 6.91, height = 2.66)
print(ggarrange(PMGraph1, PMGraph2, PMGraph3, labels = c("A","B","C"), 
          ncol = 3))
grDevices::dev.off()

# ---- Reorder & Rename --------------------------------------------------------
dnam <- dnam[, c((ncol(dnam) - 2):ncol(dnam), 1:(ncol(dnam) - 3))]
dnam$EarID <- as.factor(dnam$EarID)

# ---- Merge with phenotypic data ----------------------------------------------
data <- left_join(pheno, dnam, by = c("EarID", "Date"))
# Keep intermediate objects available for inspection; avoid clearing the workspace.
# The assembled ecological dataset is updated below with derived field variables.

# ---- Reformat columns in the dataset -----------------------------------------
names(data)[names(data) == "Sex.x"] <- "Sex" # rename sex column
data$Age <- factor(data$Age, levels = c("J","SA","A"))
data$AnthelminthicTmt <- factor(data$AnthelminthicTmt, levels = c("Control","Anthelminthic"))
data$Age2 <- plyr::revalue(data$Age, c("A" = "Adult", "SA" = "Adult", "J" = "Nonadult"))
data$Repro <- factor(data$Repro, levels = c("A","D","S","N","P","G","L"))
data$Repro2 <- factor(dplyr::recode(as.character(data$Repro),
                            A = "Non-Reproductive", D = "Reproductive", S = "Reproductive",
                            N = "Non-Reproductive", P = "Non-Reproductive", G = "Reproductive", L = "Reproductive"),
                     levels = c("Non-Reproductive", "Reproductive"))
data$Mites <- as.character(data$Mites); data$Mites <- gsub("10\\+", "10", data$Mites); data$Mites <- as.numeric(data$Mites)
data$Ecto <- data$Ticks + data$Fleas + data$Mites
data$Jul <- as.integer(format(data$Date, "%j"))
data$Jul2 <- data$Jul/365
data$Batch <- droplevels(data$Batch)
data$Year <- as.factor(format(data$Date, "%Y"))
data$ComboTmt <- factor(data$ComboTmt, levels = c("Control-Control","Control-Anthelminthic",
                                                    "Supplement-Control","Supplement-Anthelminthic"))


# IVNum: number of IV treatments a mouse received including today
#detach("package:plyr", unload = TRUE)
data <- data %>%
  dplyr::mutate(row_id = row_number()) %>%  # to preserve original order
  dplyr::mutate(IVNum = map_int(row_id, function(i) {
    this_mouse <- AnalysisID[i]
    this_date <- Date[i]
    sum(data$AnalysisID == this_mouse & data$Date <= this_date & data$TmtGiven == "IV")
  })) %>%
  dplyr::arrange(row_id) %>%
  dplyr::select(-row_id)

# TmtNum: number of C/IV treatments a mouse received including today
data <- data %>%
  dplyr::mutate(row_id = row_number()) %>%  # to preserve original order
  dplyr::mutate(TmtNum = map_int(row_id, function(i) {
    this_mouse <- AnalysisID[i]
    this_date <- Date[i]
    sum(data$AnalysisID == this_mouse & data$Date <= this_date & data$TmtGiven %in% c("C", "IV"))
  })) %>%
  dplyr::arrange(row_id) %>%
  dplyr::select(-row_id)

data$FirstCap <- ifelse(data$TmtNum == 1, "Pre-Treatment", "Post-Treatment")
data$FirstCap <- factor(data$FirstCap, levels = c("Pre-Treatment", "Post-Treatment"))

# Number of days since last IV treatment (IV mice only)
df_iv <- data %>%
  filter(TmtGiven == "IV") %>%
  arrange(AnalysisID, Date) %>%
  group_by(AnalysisID) %>%
  dplyr::mutate(
    days_since_last_IV = as.numeric(Date - lag(Date))
  ) %>%
  ungroup()

data <- data %>%
  left_join(
    df_iv %>% dplyr::select(AnalysisID, Date, days_since_last_IV),
    by = c("AnalysisID", "Date")
  )
rm(df_iv)

data <- data %>% ungroup()

# Number of days since first DNA sample (0 for first ear punch)
first_dates <- data %>%
  filter(!is.na(predWMAge)) %>%
  group_by(AnalysisID) %>%
  dplyr::summarize(first_date = min(Date), .groups = "drop")

data <- data %>%
  left_join(first_dates, by = "AnalysisID") %>%
  dplyr::mutate(
    days_since_first = if_else(
      !is.na(predWMAge),                          # only when X has a value
      as.numeric(Date - first_date),      # do the subtraction
      NA_real_                            # otherwise, NA
    )
  ) %>%
  dplyr::select(-first_date)
rm(first_dates)

# Timepoint 1 or 2 for DNAm analyses (first derivation; subsequently overwritten)
data <- data %>%
  dplyr::mutate(
    Timepoint = if_else(days_since_first == 0, "First",
                        if_else(days_since_first > 0, "Second", NA_character_))
  )

# Final Timepoint assignment below uses row order within each animal (all captures),
# not strictly the first/second DNAm observation; retained to preserve model inputs.
data <- data %>%
  dplyr::arrange(AnalysisID, Date) %>%                # make sure rows are ordered by person and time
  group_by(AnalysisID) %>%                     # group by individual
  dplyr::mutate(
    Timepoint = case_when(
      row_number() == 1 ~ "First",     # first row within each person
      row_number() == 2 ~ "Second",    # second row within each person
      TRUE ~ NA_character_             # if someone has more than 2 rows, label others NA
    ))


condition_by_year <- function(df, var1, var2, id, date_col, year_col) {
  require(lmodel2)
  require(dplyr)
  
  # Rename key variables for internal use
  df <- df %>%
    dplyr::mutate(
      v1 = !!sym(var1),
      v2 = !!sym(var2),
      MouseID = !!sym(id),
      Date = !!sym(date_col),
      Year = !!sym(year_col)
    )
  
  
  # Process each year independently
  cond_by_year <- df %>%
    group_by(Year) %>%
    group_split() %>%
    purrr::map_dfr(function(data_year) {
      # Keep original subset for later merging
      original_data <- data_year
      
      # Filter to rows with complete data for SMA and regression
      data_year <- data_year %>% filter(!is.na(v1), !is.na(v2))
      if (nrow(data_year) < 3) {
        original_data$smi <- NA_real_
        original_data$rindex <- NA_real_
        return(original_data)
      }
      
      # Compute SMA slope
      data_year <- data_year %>%
        dplyr::mutate(
          logv1 = log(v1),
          logv2 = log(v2)
        )
      
      sma <- lmodel2(logv1 ~ logv2, data = data_year, "relative", "relative", 99)
      bSMA <- sma$regression.results$Slope[sma$regression.results$Method == "SMA"]
      
      # Compute average length for standardization
      L0.avg <- data_year %>%
        group_by(MouseID) %>%
        dplyr::summarize(length = mean(v2), .groups = "drop") %>%
        dplyr::summarize(L0 = mean(length)) %>%
        pull(L0)
      
      # Calculate SMI and residual index on complete cases
      data_year <- data_year %>%
        mutate(
          smi = v1 * (L0.avg / v2)^bSMA
        )
      
      rindex_model <- lm(v1 ~ v2, data = data_year)
      data_year$rindex <- resid(rindex_model)
      
      # Join back to original group
      original_data <- original_data %>%
        left_join(data_year %>% dplyr::select(MouseID, Date, smi, rindex),
                  by = c("MouseID", "Date"))
      
      return(original_data)
    })
  
  return(cond_by_year)
}

data <- condition_by_year(data, var1 = "BodyMass", var2 = "BodyLength", 
                          id = "AnalysisID", date_col = "Date", year_col = "Year")


# Make HPoly an integer and presence/absence
data$HPolyPA <- as.factor(ifelse(data$GramHPoly > 0, "P","A"))
data$HPolyInt <- round(log(data$GramHPoly+1))

# Calculate max time observed per ID and join back
data <- data %>%
  group_by(AnalysisID) %>%
  dplyr::mutate(max_time_observed = as.numeric(max(Date) - min(Date))) %>%
  ungroup()

# Save the assembled field + DNAm dataset after calculating ecological covariates.
# This preserves all records, including samples without DNAm measurements.
save(data, file = file.path(processed_dir, "WMClock_EcologicalData.Rdata"))

################################################################################
## Samples Size Information ####################################################
################################################################################
summary_stats <- data %>%
  group_by(Year, ComboTmt) %>%
  dplyr::summarize(
    total_mice = n_distinct(AnalysisID),
    female_mice = n_distinct(AnalysisID[Sex == "F"]),
    male_mice = n_distinct(AnalysisID[Sex == "M"]),
    .groups = "drop"
  )
print(summary_stats)

################################################################################
## Does anthelminthic treatment reduce H. polygyrus egg burden? ###############
################################################################################

# Total number of faecal samples (i.e., total rows)
n_total <- nrow(subset(data,!is.na(GramHPoly)))

# Summary stats (number of samples per mouse, then summarize across mice)
fecal_counts <- data %>% subset(!is.na(HPolyInt)) %>%
  count(AnalysisID, name = "n_samples")

summary_stats <- fecal_counts %>%
  dplyr::summarize(
    total_samples = sum(n_samples),
    mean_samples = mean(n_samples),
    median_samples = median(n_samples),
    min_samples = min(n_samples),
    max_samples = max(n_samples),
    se_samples = sd(n_samples) / sqrt(n())  # standard error
  )

# Summary stats of FEC H poly burden
HPolyStats <- data %>%
  dplyr::summarize(
    nrows = nrow(data[!is.na(data$GramHPoly),]),
    mean_FEC = mean(GramHPoly, na.rm = T),
    median_FEC = median(GramHPoly, na.rm = T),
    min_FEC = min(GramHPoly, na.rm = T),
    max_FEC = max(GramHPoly, na.rm = T),
    se_FEC = sd(GramHPoly, na.rm = T) / sqrt(nrows)
  )

# ---- Calculate percentage of individuals with h poly -------------------------
hpsummary <- data %>%
  group_by(Year, AnthelminthicTmt, FirstCap) %>%
  dplyr::summarise(
    total = n_distinct(AnalysisID),
    with_P = n_distinct(AnalysisID[HPolyPA == "P"]),
    percent_P = round((with_P / total) * 100, 1),
    .groups = "drop"
  )

# ---- Anthelminthic treatment and H. polygyrus egg burden ---------------------
# Scientific question:
#   Does the difference in H. polygyrus faecal egg burden between anthelminthic
#   and control mice change between the pre- and post-treatment sampling periods?
#   The AnthelminthicTmt:FirstCap interaction is the primary treatment test.
#
# Response used for this manuscript analysis:
#   HPolyInt = round(log(GramHPoly + 1)). This is an integer-valued transformed
#   egg-per-gram measurement, NOT a direct count of eggs. The ZIP likelihood is
#   therefore a pragmatic description of this response; its distributional
#   assumptions should be acknowledged, and raw-count/EPG sensitivity analyses
#   are advisable before treating this as a definitive count-based analysis.
#
# Model choice:
#   The earlier zero-inflated negative-binomial model had a non-positive-definite
#   Hessian (NB dispersion approached the Poisson limit). A Poisson model fitted
#   successfully. Its mouse random-intercept variance was effectively zero;
#   dropping that term improved AIC by 2 without changing likelihood appreciably.
#   The no-random-effect specification reproduces the converged ZIP model.
#   Repeated observations within AnalysisID remain a limitation when estimating
#   standard errors; the model should not be interpreted as establishing
#   independence of repeated observations.
#
# IMPORTANT: na.omit (rather than na.exclude) is needed here. With na.exclude,
#   DHARMa returned 827 simulated/observed values but 1245 predictions, causing
#   the error 'variable lengths differ (found for pred)' when plotting.
WormMod <- glmmTMB::glmmTMB(
  HPolyInt ~ SupplementTmt + AnthelminthicTmt * FirstCap +
    Jul2 + Repro2 + Sex + Site + Age2 + Year,
  ziformula = ~ AnthelminthicTmt,
  family = poisson(link = "log"),
  data = data,
  na.action = na.omit
)

print(summary(WormMod))
print(car::Anova(WormMod, type = 2))

# Do not use a fitted model with a failed Hessian or optimizer status.
if (WormMod$fit$convergence != 0 || !isTRUE(WormMod$sdr$pdHess)) {
  stop("H. polygyrus ZIP model did not converge reliably; inspect diagnose(WormMod).")
}

# Simulation-based model diagnostics; displayed in RStudio, not saved as a
# manuscript figure. A nonsignificant test does not prove model adequacy.
set.seed(123)
worm_dharma <- DHARMa::simulateResiduals(fittedModel = WormMod, n = 500)
stopifnot(
  length(worm_dharma$observedResponse) == nobs(WormMod),
  length(worm_dharma$fittedPredictedResponse) == nobs(WormMod),
  length(worm_dharma$scaledResiduals) == nobs(WormMod)
)
plot(worm_dharma)
print(DHARMa::testUniformity(worm_dharma, plot = FALSE))
print(DHARMa::testDispersion(worm_dharma, plot = FALSE))
print(DHARMa::testZeroInflation(worm_dharma, plot = FALSE))
print(DHARMa::testOutliers(worm_dharma, plot = FALSE))

# Interpret the interaction using the *combined expected response* (conditional
# Poisson mean multiplied by 1 minus the structural-zero probability).
# These estimated means are on the HPolyInt scale, NOT original eggs per gram.
worm_emm <- emmeans::emmeans(
  WormMod,
  ~ AnthelminthicTmt * FirstCap,
  component = "response"
)
print(summary(worm_emm, infer = TRUE))

# Contrasts of treated vs control mice separately before/after treatment.
print(pairs(worm_emm, by = "FirstCap"))

# Primary marginal difference-in-differences contrast. Negative estimate:
# the treated-vs-control difference becomes more negative after treatment.
worm_interaction <- emmeans::contrast(
  worm_emm,
  interaction = c("revpairwise", "revpairwise")
)
print(summary(worm_interaction, infer = TRUE))

# ---- Graph of IV treatment effect on h poly numbers --------------------------
hpgraph <- ggplot(subset(data, !is.na(GramHPoly)), aes(x=SupplementTmt, y=log(GramHPoly+1), fill=AnthelminthicTmt)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(),aes(group=AnthelminthicTmt), col="black", pch=1) + 
  theme_classic() +
  labs(x="Food Supplementation", y = "Faecal Egg Burden (log(x+1) eggs/gram)") + 
  scale_fill_manual(name = "Anthelminthic Treatment", values = c("gray94","coral")) + 
  theme(text=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title = element_text(size=12),
        legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10)) + 
  facet_wrap(~FirstCap)
hpgraph

# save graph as a PDF
grDevices::pdf(file.path(figures_dir, "Figure_S5.pdf"), width = 5.05, height = 4.61)
print(hpgraph)
grDevices::dev.off()

################################################################################
## Does smi (skeletal muscle) and rindex change with treatment? ################
################################################################################

# ---- Residual Index ----------------------------------------------------------
hist(data$rindex)
mbc <- lmer(rindex ~ 
              SupplementTmt + AnthelminthicTmt +
              Jul2 + Repro2 + Sex + Site + Age2 + Year + 
              (1|AnalysisID),
            data = data, na.action=na.exclude)
summary(mbc)
car::Anova(mbc, type = 2, test.statistic = "F")
hist(resid(mbc))

mbc <- ggplot(data, aes(x=AnthelminthicTmt, y=rindex, fill=SupplementTmt)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(),aes(group=SupplementTmt), col="black", pch=1) + 
  theme_classic() +
  labs(x="Anthelminthic Treatment", y = "rindex")+
  scale_fill_brewer(palette="RdBu", direction =-1) +
  theme(text=element_text(size=12), legend.text=element_text(size=12), 
        axis.text=element_text(size=12), axis.title = element_text(size=12), 
        legend.title=element_blank(), legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10))
mbc

# ---- Skeletal Muscle Index ---------------------------------------------------
hist(data$smi)
sbc <- lmer(smi ~ 
              SupplementTmt + AnthelminthicTmt +
              Jul2 + Repro2 + Sex + Site + Age2 + Year +
              (1|AnalysisID),
            data = data, na.action=na.exclude)
summary(sbc)
car::Anova(sbc, type = 2, test.statistic = "F")
hist(resid(sbc))

sbc <- ggplot(data, aes(x=AnthelminthicTmt, y=smi, fill=SupplementTmt)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(),aes(group=SupplementTmt), col="black", pch=1) + 
  theme_classic() +
  labs(x="Anthelminthic Treatment", y = "skeletal muscle index")+
  scale_fill_brewer(palette="RdBu", direction =-1) +
  theme(text=element_text(size=12), legend.text=element_text(size=12), 
        axis.text=element_text(size=12), axis.title = element_text(size=12), 
        legend.title=element_blank(), legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10))
sbc

# ---- Body Mass ---------------------------------------------------------------
hist(data$BodyMass)
bbc <- lmer(BodyMass ~ 
              SupplementTmt + AnthelminthicTmt +
              Jul2 + Repro2 + Sex + Site + Age2 + Year +
              (1|AnalysisID),
            data = data, na.action=na.exclude)
summary(bbc)
car::Anova(bbc, type = 2, test.statistic = "F")
hist(resid(bbc))

bbc <- ggplot(data, aes(x=AnthelminthicTmt, y=BodyMass, fill=SupplementTmt)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(),aes(group=SupplementTmt), col="black", pch=1) + 
  theme_classic() +
  labs(x="Anthelminthic Treatment", y = "Body Mass (g)")+
  scale_fill_brewer(palette="RdBu", direction =-1) +
  theme(text=element_text(size=12), legend.text=element_text(size=12), 
        axis.text=element_text(size=12), axis.title = element_text(size=12), 
        legend.title=element_blank(), legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10))
bbc


################################################################################
## Does reproduction change with treatment? ####################################
################################################################################
rep <- glmer(Repro2 ~ 
              SupplementTmt + AnthelminthicTmt +
              SupplementTmt * Sex +
              Jul2 + Sex + Site + Year +
              (1|AnalysisID),
              data = data[!is.na(data$Sex) & !is.na(data$Repro2),], family = "binomial", 
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)), 
              na.action=na.exclude)
summary(rep)
simulateResiduals(rep) %>% plot()

################################################################################
## Does MNKA change with treatment? ############################################
################################################################################

# ---- Calculate the MNKA per Year, Treatment, and Session ---------------------
# Unique grid-year combinations
grid_years <- unique(data[, c("Grid", "Year")])

# Updated MNA calculation function for a subset of data
calc_mna <- function(data_sub) {
  unique_dates <- unique(data_sub$Session)
  
  mna_data <- sapply(unique_dates, function(x) {
    current_date <- x
    caught_at_date <- unique(data_sub$AnalysisID[data_sub$Session == current_date])
    caught_aft <- unique(data_sub$AnalysisID[data_sub$Session > current_date])
    caught_bf <- unique(data_sub$AnalysisID[data_sub$Session < current_date])
    caught_bf_and_aft <- setdiff(intersect(caught_bf, caught_aft), caught_at_date)
    length(caught_at_date) + length(caught_bf_and_aft)
  })
  
  data.frame(Session = unique_dates, Count = mna_data)
}

# Loop through each grid-year pair
long_results <- do.call(rbind, apply(grid_years, 1, function(row) {
  grid <- row[["Grid"]]
  year <- row[["Year"]]
  
  data_sub <- data %>% filter(Grid == grid, Year == year)
  
  if (nrow(data_sub) == 0) return(NULL)  # skip empty combos
  
  combo_tmt <- unique(data_sub$SupplementTmt)
  if (length(combo_tmt) != 1) {
    warning(paste("Grid-Year has multiple treatments:", grid, year))
    combo_tmt <- combo_tmt[1]  # handle ambiguity by selecting first
  }
  
  mna_df <- calc_mna(data_sub)
  mna_df$Grid <- grid
  mna_df$Year <- year
  mna_df$Treatment <- combo_tmt
  return(mna_df)
}))

long_results$Site <- factor(dplyr::recode(as.character(long_results$Grid),"1" = "Hewan", "2" = "Hewan", 
                     "7" = "Penicuik", "9" = "Penicuik"), levels = c("Penicuik","Hewan"))
                             
# ---- Prepare a graph ---------------------------------------------------------
MNKAgraph <- ggplot(subset(long_results, !is.na(Site)), aes(x = Session, y = Count, colour = Treatment, group = Grid)) + 
  geom_point() + 
  geom_line(aes(lty = Site)) + 
  facet_wrap(~Year) + 
  theme_bw() + 
  labs(x = "Trapping Session", y = "Min # Known Alive (per grid)") + 
  scale_color_manual(values = c("black","orange")) + 
  scale_shape_manual(values = c(16,1)) + 
  scale_x_continuous(limits = c(01,9), breaks = 1:9)
MNKAgraph

# save graph as a PDF
grDevices::pdf(file.path(figures_dir, "Figure_S4.pdf"), width = 5.14, height = 3.59)
print(MNKAgraph)
grDevices::dev.off()

# ---- Test whether MNKA differs by food treatment -----------------------------
MNKA <- glm(Count ~ Session*Treatment + Year + Site, data=long_results, family = "poisson")
summary(MNKA)            
simulateResiduals(MNKA) %>% plot()

# --- See whether the number of juveniles seen in each treatment differ --------
# Count unique juveniles per ComboTmt and Year
juvenile_counts <- data %>%
  subset(Age == "J") %>%
  group_by(Year, SupplementTmt) %>%
  summarize(Unique_Juveniles = n_distinct(AnalysisID), .groups = "drop")

################################################################################
## Does epigenetic age fit with age categories? ################################
################################################################################
# All time points per mouse
agecat <- lmer(predWMAge ~ Age + (1|AnalysisID),
                data = data, na.action = na.omit)
summary(agecat)
car::Anova(agecat, type = 2, test.statistic = "F")
nobs(agecat)
hist(resid(agecat))

# posthoc test of age category
emm <- emmeans(agecat, ~ Age)
contrast(emm, method = "revpairwise")

# age categories
AgeCat <- ggplot(subset(data, !is.na(Age)), aes(x=Age,y=predWMAge)) + 
  #geom_point(data=rando, aes(x=Cat, y=predAgeAll), color = "red",pch=2)+
  geom_boxplot() + 
  geom_point(aes(x=Age, y=predWMAge, fill=Age),position=position_jitterdodge(jitter.width = 0.75),
             pch = 1, size = 2, color = "black") + 
  theme_bw() + 
  labs(x="Morphological Age Category", y="") + 
  scale_x_discrete(labels=c("Juvenile","Subadult","Adult")) + 
  theme(legend.position = "none")
AgeCat

# summaries
catsum <- data %>%
  filter(!is.na(Age), !is.na(predWMAge)) %>%
  group_by(Age) %>%
  dplyr::summarize(
    nrows = n(),
    meane = mean(predWMAge, na.rm = T),
    mediane = median(predWMAge, na.rm = T),
    mine = min(predWMAge, na.rm = T),
    maxe = max(predWMAge, na.rm = T),
    see = sd(predWMAge, na.rm = T) / sqrt(nrows)
  )
print(catsum)

catsum2 <- data %>%
  filter(!is.na(Age), !is.na(predWMAge)) %>%
  dplyr::summarize(
    nrows = n(),
    meane = mean(predWMAge, na.rm = T),
    mediane = median(predWMAge, na.rm = T),
    mine = min(predWMAge, na.rm = T),
    maxe = max(predWMAge, na.rm = T),
    see = sd(predWMAge, na.rm = T) / sqrt(nrows)
  )
print(catsum2)

# save graph as a PDF
# Plot displayed in RStudio but not exported: AgeCat.pdf

# age categories juveniles vs non-juveniles
AgeCatJ <- ggplot(subset(data, !is.na(Age2)), aes(x=Age2,y=predWMAge)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(x=Age2, y=predWMAge, fill=Age2),position=position_jitterdodge(jitter.width = 0.75),
             pch = 1, size = 2, color = "black") + 
  theme_bw() + 
  labs(x="Morphological Age Category", y="Predicted Age (Days)") + 
  scale_x_discrete(labels=c("Juvenile","Non-Juvenile")) + 
  theme(legend.position = "none")
AgeCatJ

################################################################################
## What predicts epigenetic age? ###############################################
################################################################################

# only first time point per mouse
temp <- data %>% subset(FirstCap == "Pre-Treatment") %>%
  subset(!is.na(predWMAge))
epiAge <- lm(predWMAge ~ SupplementTmt + Site + Jul2 + 
                 Sex + Repro2 + BodyMass + Year + Age2,
               data = temp, na.action = na.omit)
summary(epiAge)
car::vif(epiAge)
car::Anova(epiAge, type = 2, test.statistic = "F")
visreg(epiAge, "SupplementTmt")
nobs(epiAge)

# All time points per mouse (REPORT THIS ONE)
epiAge2 <- lmer(predWMAge ~ SupplementTmt + AnthelminthicTmt + Site + Jul2 + Timepoint + 
               Sex + Repro2 + BodyMass + Year + 
                 (1|AnalysisID),
             data = data, na.action = na.omit)
summary(epiAge2)
car::vif(epiAge2)
car::Anova(epiAge2, type = 2, test.statistic = "F")
nobs(epiAge2)

massgraph <- ggplot(data, aes(x = BodyMass, y = predWMAge)) + 
  geom_point(pch = 1, size = 2) + 
  geom_smooth(method = "lm", se = F, color = "black") + 
  labs(x = 'Body Mass (grams)', y = '') + 
  theme_bw()
massgraph

data$DateConvert <- as.Date(data$Jul - 1, origin = "2024-01-01")
julgraph <- ggplot(data, aes(x = DateConvert, y = predWMAge)) + 
  geom_smooth(aes(group = AnalysisID), method = "lm", se = F, color = "gray50") + 
  geom_point(pch = 16, size = 2) + 
  labs(x = 'Date', y = '') + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month")  + 
  theme_bw()
julgraph

julgraph2 <- ggplot(data, aes(x = DateConvert, y = predWMAge, fill = AnalysisID)) + 
  geom_smooth(method = "lm", se = F, color = "gray") + 
  geom_point(pch = 16, size = 2) + 
  labs(x = 'Date', y = '') + 
  scale_x_date(date_labels = "%b", date_breaks = "1 month")  + 
  theme_bw() + 
  theme(legend.position = "none")
julgraph2

# Make sure Julian0 exists in your actual data
data$DateConvert <- as.Date(data$Jul - 1, origin = "2024-01-01")
# Compute Julian0
data <- data %>%
  mutate(Julian0 = as.numeric(DateConvert - min(DateConvert)))

# Month breaks and labels
month_breaks <- seq(
  from = floor_date(min(data$DateConvert), "month"),
  to = ceiling_date(max(data$DateConvert), "month"),
  by = "1 month"
)
month_breaks_julian <- as.numeric(month_breaks - min(data$DateConvert))
month_labels <- format(month_breaks, "%b")

# Remove the last month (May) from labels
month_breaks_julian_plot <- month_breaks_julian[-length(month_breaks_julian)]
month_labels_plot <- month_labels[-length(month_labels)]

# Bottom positions for ticks and labels (safe for 0-400 y-range)
tick_y <- 0 + 0.01 * 400
label_y <- 0 + 0.03 * 400

slopes <- data %>%
  group_by(AnalysisID) %>%
  filter(!is.na(predWMAge), !is.na(Julian0)) %>%    # remove NA rows
  filter(n() >= 2) %>%                              # at least 2 points to fit a line
  dplyr::summarise(
    slope = coef(lm(predWMAge ~ Julian0))[2],
    .groups = "drop"
  ) %>%
  dplyr::mutate(slope_sign = ifelse(slope >= 0, "positive", "negative"))

data <- data %>%
  dplyr::select(-any_of("slope_sign")) %>%  # allows this section to be rerun safely
  left_join(slopes %>% select(AnalysisID, slope_sign), by = "AnalysisID")

# Plot
# Plot
julgraph3 <- ggplot(data, aes(x = Julian0, y = predWMAge)) +
  
  # Within-individual trajectories
  geom_smooth(
    data = data %>% filter(!is.na(slope_sign)),
    aes(color = slope_sign, group = AnalysisID),
    method = "lm",
    se = FALSE,
    linewidth = 0.8
  ) +
  
  # Individual observations
  geom_point(
    pch = 16,
    size = 2
  ) +
  
  # Population-level seasonal trend
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 2
  ) +
  
  scale_color_manual(
    values = c(
      "positive" = "#FF6666",
      "negative" = "#6699FF"
    )
  ) +
  
  scale_x_continuous(
    name = "Days Since Season Start",
    limits = c(0, max(data$Julian0, na.rm = TRUE))
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(y = "") +
  
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  
  annotate(
    "text",
    x = month_breaks_julian_plot,
    y = label_y,
    label = month_labels_plot,
    vjust = 0.75,
    size = 3.5
  )

julgraph3

# save graph as a PDF
# Plot displayed in RStudio but not exported: UpdatedMonthGraph.pdf

foodgraph <- ggplot(data, aes(x=SupplementTmt, y=predWMAge, fill=AnthelminthicTmt)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(),aes(group=AnthelminthicTmt), col="black", pch=1, size = 2) + 
  theme_bw() +
  labs(x="Food Supplementation", y = "Predicted Epigenetic Age") + 
  scale_fill_manual(name = "", values = c("gray94","coral")) + 
  theme(text=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title = element_text(size=12),
        legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10))
foodgraph

reprograph <- ggplot(data[!is.na(data$Repro2) & data$Sex=="F",], aes(x=Repro2, y=predWMAge)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(),
             aes(color = BodyMass), pch=16, size = 2) + 
  theme_bw() +
  labs(x="Reproductive Status", y = "Predicted Epigenetic Age") + 
  scale_color_gradient2(low = "red", mid = "green", high = "blue",
                        midpoint = mean(data$BodyMass, na.rm = TRUE)) + 
  theme(text=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title = element_text(size=12),
        legend.position = "bottom",
        legend.position.inside = c(0.3,0.1),
        plot.margin = ggplot2::margin(10,10,10,10))
reprograph

# save graph as a PDF
# Plot displayed in RStudio but not exported: EpiAgeCorrelates.pdf

################################################################################
# Is epiAge related to parasites? ##############################################
################################################################################
# reformat blood volume column to be numerical
data$BloodNum <- as.character(data$BloodVol)
data$BloodNum[data$BloodNum == ">50"]<-"50"
data$BloodNum[data$BloodNum == "<1"]<-"1"
data$BloodNum <- as.numeric(data$BloodNum)

# ---- Does EpiAge predict the number of types of parasites you have? ----------
# Create new Y/N columns for each type and then a column with their sum of Y's
data <- data %>%
  dplyr::mutate(
    HPolyY    = GramHPoly > 0,
    TicksY    = Ticks > 0,
    FleaseY   = Fleas > 0,
    MitesY    = Mites > 0,
    WMHVY     = WMHV.PA == "P",
    trypsY    = Tryps.PA == "P",
    bartsY    = Barts.PA == "P"
  ) %>%
  dplyr::mutate(across(c(HPolyY, TicksY, FleaseY, MitesY, WMHVY, trypsY, bartsY), as.integer)) %>%
  dplyr::mutate(ParCount = rowSums(across(c(HPolyY, TicksY, FleaseY, MitesY, WMHVY, trypsY, bartsY))))


# Tried this with the random effect, but it explained 0% of the variance, so removed.
# Model - all timepoints per mouse
# PNum <- glm(ParCount ~ SupplementTmt + AnthelminthicTmt + predWMAge + Jul2 + 
#             Sex + Repro2 + Age2 + Year, 
#             data = data, family = "poisson",
#             na.action = na.omit)
# summary(PNum)
# car::vif(PNum)
# simulateResiduals(PNum) %>% plot()

#Model - only first timepoint per mouse
temp <- data %>% subset(Timepoint == "First")
PNum <- glm(ParCount ~ SupplementTmt + predWMAge + Jul2 + Sex + Repro2 + Year,
                data = temp, family = "poisson",
                na.action = na.omit)
summary(PNum)
nobs(PNum)
car::vif(PNum)
simulateResiduals(PNum) %>% plot()

# Create a graph
ParNum <- ggplot(temp, aes(x=ParCount, y=predWMAge)) + 
  geom_point(pch=1, size = 2.5, position = position_jitter(width = 0.1)) +
  geom_smooth(method="lm", se = T, color = "black") + 
  theme_bw() + 
  lims(x = c(0,6)) + 
  labs(x = "Parasite Species Richness", y = "Predicted Age (Days)") + 
  theme(legend.position = "none")
ParNum

# save graph as a PDF
# Plot displayed in RStudio but not exported: ParasiteTypeCount.pdf

# ---- Is epigenetic age predictive of specific parasites? ---------------------
Pre <- data %>% subset(FirstCap == "Pre-Treatment")
Post <- data %>% subset(FirstCap == "Post-Treatment")

WMHV_plot <- ggplot(Pre[!is.na(Pre$WMHV.PA),], aes(x = WMHV.PA, y = predWMAge)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(x = WMHV.PA, y = predWMAge, fill = WMHV.PA), 
  alpha = 0.6, position = position_jitterdodge(jitter.width = 0.5),
  pch = 1, size = 2.5, color = "black") + 
  theme_bw() +
  labs(x = "Wood Mouse Herpes Virus", y = "") +
  scale_x_discrete(labels = c("A" = "Absent", "P" = "Present")) +
  theme(legend.position = "none")

Bart_plot <- ggplot(Pre[!is.na(Pre$Barts.PA),], aes(x = Barts.PA, y = predWMAge)) +
  geom_boxplot() + 
  geom_point(aes(x = Barts.PA, y = predWMAge, fill = Barts.PA), 
             alpha = 0.6, position = position_jitterdodge(jitter.width = 0.5),
             pch = 1, size = 2, color = "black") + 
  theme_bw() +
  labs(x = "Bartonella Presence", y = "WM EpiAge") +
  theme(legend.position = "none")

Tryps_plot <- ggplot(Pre[!is.na(Pre$Tryps.PA),], aes(x = Tryps.PA, y = predWMAge)) +
  geom_boxplot() + 
  geom_point(aes(x = Tryps.PA, y = predWMAge, fill = Tryps.PA), 
             alpha = 0.6, position = position_jitterdodge(jitter.width = 0.5),
             pch = 1, size = 2, color = "black") + 
  theme_bw() +
  labs(x = "Trypanosoma Presence", y = "WM EpiAge") +
  theme(legend.position = "none")

Ticks_plot <- ggplot(Pre, aes(x = log(Ticks + 1), y = predWMAge)) +
  geom_point(alpha = 0.6, pch = 1, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  labs(x = "Log(Number of Ticks + 1)", y = "") +
  theme(legend.position = "right")

Mites_plot <- ggplot(Pre, aes(x = log(Mites), y = predWMAge)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(x = "# Mites", y = "WM EpiAge") +
  theme(legend.position = "right")

Fleas_plot <- ggplot(Pre, aes(x = log(Fleas), y = predWMAge)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(x = "# Fleas", y = "WM EpiAge") +
  theme(legend.position = "right")

Hpoly_plot <- ggplot(Pre, aes(x = HPolyInt, y = predWMAge)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(x = "log(# HPoly)", y = "WM EpiAge") +
  theme(legend.position = "right")

# Put all plots together into one
ggarrange(Ticks_plot, Mites_plot, Fleas_plot, Hpoly_plot, WMHV_plot, Bart_plot, 
          Tryps_plot, nrow=2,ncol=4)


# ---- Multivariate analyses of parasites --------------------------------------
multi <- data[,c("predWMAge","SupplementTmt","Jul2","Sex","Repro2","Age2","Year",
                 "Ticks","Mites","Fleas","WMHV.PA","Tryps.PA",
                 "Barts.PA","HPolyInt","Mouse.ID","Timepoint")]
multi <- multi[complete.cases(multi), ]
multi <- multi %>% subset(Timepoint == "First")
multi <- as.data.frame(multi)
multi <- droplevels(multi)


manova_model <- manova(cbind(Ticks, Mites, Fleas, Barts.PA, WMHV.PA, Tryps.PA, HPolyInt) 
                       ~ predWMAge + SupplementTmt + Jul2 + Sex + Repro2 + Year,
                       data = multi)
summary(manova_model)
# Summary of univariate ANOVAs
mansum <- summary.aov(manova_model)

fvals <- sapply(mansum, function(tbl) tbl[["F value"]])
pvals <- sapply(mansum, function(tbl) tbl[["Pr(>F)"]])

all_p <- unlist(pvals)
all_p_adj <- p.adjust(all_p, method = "fdr")

# Convert adjusted p-values back into a matrix
pvals_adj <- matrix(all_p_adj, nrow = nrow(pvals), byrow = FALSE)
rownames(pvals_adj) <- rownames(pvals)
colnames(pvals_adj) <- colnames(pvals)

# Extract F and p-values into a tidy format
options(scipen = 999)
results_list <- lapply(names(mansum), function(resp) {
  df <- mansum[[resp]]
  data.frame(
    Response = resp,
    Predictor = rownames(df),
    F_value = df[,"F value"],
    p_value = df[,"Pr(>F)"]
  )
})

# Combine all response variables
big_results <- bind_rows(results_list)

# Adjust p-values across all tests
big_results$p_adj <- p.adjust(big_results$p_value, method = "fdr")
big_results <- big_results %>%
  mutate(
    F_value = round(F_value, 4),
    p_value = round(p_value, 4),
    p_adj = round(p_adj, 4)
  )
big_results

# save all parasite-related graphs as a PDF
# Plot displayed in RStudio but not exported: RelevantParasiteGraphs.pdf

################################################################################
## Longitudinal Analyses of Epigenetic Age #####################################
################################################################################
# Package conflict is handled with explicit dplyr:: calls where required.

WildLong <- data %>%
  subset(!is.na(predWMAge)) %>%
  group_by(AnalysisID) %>%
  summarize(
    first_epi = predWMAge[which.min(Date)],
    last_epi = predWMAge[which.max(Date)],
    first_date = min(Date),
    last_date = max(Date),
    Food = SupplementTmt[which.min(Date)],
    IV = AnthelminthicTmt[which.min(Date)],
    Combo = ComboTmt[which.min(Date)],
    TmtNum = TmtNum[which.max(Date)],
    Sex = Sex[which.min(Date)],
    Age = Age[which.min(Date)],
    Repro.1 = Repro2[which.min(Date)],
    Repro.2 = Repro2[which.max(Date)],
    Year = Year[which.min(Date)],
    TObs = max_time_observed[which.min(Date)]
  )

reprotime <- data %>%
  # keep only rows where Repro2 has a value
  filter(!is.na(Repro2) & Repro2 != "") %>%
  group_by(AnalysisID) %>%
  summarise(
    total_sessions = n(),
    repro_sessions = sum(Repro2 == "Reproductive", na.rm = TRUE),
    percent_matching = 100 * repro_sessions / total_sessions,
    .groups = "drop"
  )
print(reprotime)

WildLong <- left_join(WildLong, reprotime, by = "AnalysisID")

WildLong$Diff <- WildLong$last_epi - WildLong$first_epi
WildLong$Dur <- as.numeric(WildLong$last_date - WildLong$first_date)
WildLong$Slope <- (WildLong$last_epi - WildLong$first_epi)/as.numeric(WildLong$last_date - WildLong$first_date)
hist(WildLong$Diff)
hist(WildLong$Slope)

WildLong$Combo <- factor(WildLong$Combo, levels = c("Control-Control","Control-Anthelminthic",
                                                    "Supplement-Control","Supplement-Anthelminthic"))

## two individuals have weird negative slopes less than -1
## 041F000599 looks like their values are backwards first vs second time point
## 041F000293 decreases by 223 in 80ish days, something seems wrong
## WildLong <- WildLong[!(WildLong$AnalysisID %in% c("041F000293", "041F0002DE")), ]

WildLong$Direction <- ifelse(WildLong$Slope > 0, "Up","Down")
WildLong$ReproChange <- as.factor(paste0(WildLong$Repro.1,"-",WildLong$Repro.2))

summary_stats <- WildLong %>%
  subset(!is.na(Slope)) %>%
  group_by(Year, Combo) %>%
  summarize(
    total_mice = n_distinct(AnalysisID),
    female_mice = n_distinct(AnalysisID[Sex == "F"]),
    male_mice = n_distinct(AnalysisID[Sex == "M"]),
    .groups = "drop"
  )
print(summary_stats)

WildLong2 <- WildLong %>% subset(Dur>35)
# WildLong: one first-to-last summary per mouse with DNAm data.
# WildLong2: longitudinal modelling subset with >35 days between samples.
save(WildLong, WildLong2,
     file = file.path(processed_dir, "WMClock_LongitudinalData.Rdata"))
stderr <- function(x) sd(x) / sqrt(length(x))
cat(mean(WildLong2$Dur), " ", stderr(WildLong2$Dur), " ", max(WildLong2$Dur), "\n")

stderr <- function(x) sd(x) / sqrt(length(x))
cat(mean(WildLong2$Diff), " ", stderr(WildLong2$Diff), " ", min(WildLong2$Diff), " ", max(WildLong2$Diff), "\n")

stderr <- function(x) sd(x) / sqrt(length(x))
cat(mean(WildLong2$Slope), " ", stderr(WildLong2$Slope), " ", min(WildLong2$Slope), " ", max(WildLong2$Slope), "\n")

sum(WildLong$Slope<0, na.rm = T)/sum(WildLong$Slope>0, na.rm = T)

## option 1. does change in day match change in epiAge?
diffmod <- lm(Diff ~ Food*IV + IV*Dur + Food*Dur + Year + Sex + Repro.1, 
              data = WildLong2, na.action=na.omit)

summary(diffmod)
#hist(resid(diffmod))
anova(diffmod)
visreg(diffmod, "IV")
nobs(diffmod)

# Get simple slopes of Dur within each IV
slopes <- emtrends(diffmod, ~ IV, var = "Dur")
summary(slopes)
# Compare those slopes across IV levels
pairs(slopes, adjust = "tukey")

# Get simple slopes of Dur within each Food
slopes <- emtrends(diffmod, ~ Food, var = "Dur")
summary(slopes)
# Compare those slopes across food levels
pairs(slopes, adjust = "tukey")

DiffplotIV <- ggplot(data = WildLong2, aes(x = Dur, y = Diff, fill = IV)) + 
  geom_point(aes(shape = IV, color = IV), size = 2, show.legend = F) +
  geom_smooth(aes(fill = IV, color = IV), method = "lm", se = T) + 
  theme_bw() + 
  labs(x = "Duration between Paired Samples", y = "Difference in EpiAge") + 
  scale_fill_manual(name = "", values = c("gray50","coral")) +
  scale_color_manual(name = "", values = c("gray50","coral")) +
  #lims(x = c(-25,400), y = c(-25,400)) + 
  geom_abline(intercept = 0, slope = 1, size = 0.75, color = "black", linetype = "dashed") + 
  theme(legend.position = "bottom")
DiffplotIV

DiffplotFood <- ggplot(data = WildLong2, aes(x = Dur, y = Diff, fill = Food)) + 
  geom_point(aes(shape = Food, color = Food), size = 2, show.legend = F) +
  geom_smooth(aes(fill = Food, color = Food), method = "lm", se = T) + 
  theme_bw() + 
  labs(x = "Duration between Paired Samples", y = "Difference in EpiAge") + 
  scale_fill_manual(name = "", values = c("gray50","red")) +
  scale_color_manual(name = "", values = c("gray50","red")) +
  #lims(x = c(-25,400), y = c(-25,400)) + 
  geom_abline(intercept = 0, slope = 1, size = 0.75, color = "black", linetype = "dashed") + 
  theme(legend.position = "bottom")
DiffplotFood

# combo treatment
diffmod <- lm(Diff ~ Combo*Dur + Year + Sex + Repro.1, 
              data = WildLong2, na.action=na.omit)

summary(diffmod)
#hist(resid(diffmod))
anova(diffmod)

# Get simple slopes of Dur
slopes <- emtrends(diffmod, ~ Combo, var = "Dur")
summary(slopes)
pairs(slopes, adjust = "tukey")

# Preserve the primary model object before running reviewer sensitivity analyses.
diffmod_primary <- diffmod
slopes_primary <- emtrends(diffmod_primary, ~ Combo, var = "Dur")

################################################################################
# Reviewer 1: sensitivity to extreme negative epigenetic-ageing slopes
################################################################################
# The original exploratory notes identified a small number of unusually large
# negative slopes. Rather than hard-coding individual IDs, use the pre-specified
# criterion already noted in this script (Slope < -1) to identify them.

extreme_negative <- WildLong2 %>%
  filter(Slope < -1) %>%
  arrange(Slope) %>%
  dplyr::select(
    AnalysisID, Combo, Year,
    first_epi, last_epi,
    Dur, Diff, Slope
  )

print(extreme_negative)
cat("Number of extreme negative slopes (< -1): ", nrow(extreme_negative), "\n", sep = "")

write.csv(
  extreme_negative,
  file.path(tables_dir, "Sensitivity_extreme_negative_individuals.csv"),
  row.names = FALSE
)

# Repeat the exact primary Combo x Duration model after excluding these animals.
WildLong_no_extreme <- WildLong2 %>%
  filter(Slope >= -1 | is.na(Slope))

cat("Primary longitudinal dataset: ", nrow(WildLong2), " mice\n", sep = "")
cat("After excluding extreme negative slopes: ",
    nrow(WildLong_no_extreme), " mice\n", sep = "")

diffmod_no_extreme <- lm(
  Diff ~ Combo * Dur + Year + Sex + Repro.1,
  data = WildLong_no_extreme,
  na.action = na.omit
)

summary(diffmod_no_extreme)
anova(diffmod_no_extreme)

slopes_no_extreme <- emtrends(
  diffmod_no_extreme,
  ~ Combo,
  var = "Dur"
)

summary(slopes_no_extreme)
pairs(slopes_no_extreme, adjust = "tukey")

# Side-by-side table of treatment-specific Dur slopes.
outlier_sensitivity_slopes <- bind_rows(
  as.data.frame(slopes_primary) %>%
    mutate(Analysis = "Primary"),
  as.data.frame(slopes_no_extreme) %>%
    mutate(Analysis = "Excluding slopes < -1")
)

print(outlier_sensitivity_slopes)
write.csv(
  outlier_sensitivity_slopes,
  file.path(tables_dir, "Sensitivity_outlier_epigenetic_ageing_slopes.csv"),
  row.names = FALSE
)

################################################################################
# Reviewer 2: sample size and elapsed-time distribution by treatment group
################################################################################
# Quantify the visual concern that the untreated control group may contain fewer
# animals and shorter intervals between paired DNAm samples.

duration_by_treatment <- WildLong2 %>%
  group_by(Combo) %>%
  summarise(
    n = n(),
    mean_Dur = mean(Dur, na.rm = TRUE),
    SE_Dur = sd(Dur, na.rm = TRUE) / sqrt(sum(!is.na(Dur))),
    median_Dur = median(Dur, na.rm = TRUE),
    min_Dur = min(Dur, na.rm = TRUE),
    Q1_Dur = as.numeric(quantile(Dur, 0.25, na.rm = TRUE)),
    Q3_Dur = as.numeric(quantile(Dur, 0.75, na.rm = TRUE)),
    max_Dur = max(Dur, na.rm = TRUE),
    n_under_50 = sum(Dur < 50, na.rm = TRUE),
    percent_under_50 = 100 * mean(Dur < 50, na.rm = TRUE),
    .groups = "drop"
  )

print(duration_by_treatment)
write.csv(
  duration_by_treatment,
  file.path(tables_dir, "Duration_by_treatment_group.csv"),
  row.names = FALSE
)

# Diagnostic plot of elapsed time by treatment group. This is intended primarily
# for checking the reviewer concern; it does not need to be a manuscript figure.
DurationByTreatment_plot <- ggplot(WildLong2, aes(x = Combo, y = Dur)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 2) +
  theme_bw() +
  labs(
    x = "Treatment",
    y = "Days between DNAm samples"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

DurationByTreatment_plot

# Plot displayed in RStudio but not exported: DurationByTreatment_diagnostic.pdf

################################################################################
# Reviewer 2: common-support sensitivity analysis for elapsed time
################################################################################
# Restrict the model to the range of Dur represented in all four treatment
# groups. This avoids choosing an arbitrary cut-off based on the reviewer's
# visual impression alone.

duration_ranges <- WildLong2 %>%
  group_by(Combo) %>%
  summarise(
    min_Dur = min(Dur, na.rm = TRUE),
    max_Dur = max(Dur, na.rm = TRUE),
    .groups = "drop"
  )

print(duration_ranges)

common_min <- max(duration_ranges$min_Dur)
common_max <- min(duration_ranges$max_Dur)

cat(
  "Common duration range across all treatment groups: ",
  common_min, " to ", common_max, " days\n",
  sep = ""
)

WildLong_common <- WildLong2 %>%
  filter(Dur >= common_min, Dur <= common_max)

common_support_n <- WildLong_common %>%
  count(Combo, name = "n")

print(common_support_n)

write.csv(duration_ranges,
          file.path(tables_dir, "Duration_range_by_treatment_group.csv"), row.names = FALSE)
write.csv(common_support_n,
          file.path(tables_dir, "Common_support_sample_sizes.csv"), row.names = FALSE)

# Repeat the same primary model in the common-support subset.
diffmod_common <- lm(
  Diff ~ Combo * Dur + Year + Sex + Repro.1,
  data = WildLong_common,
  na.action = na.omit
)

summary(diffmod_common)
anova(diffmod_common)

slopes_common <- emtrends(
  diffmod_common,
  ~ Combo,
  var = "Dur"
)

summary(slopes_common)
pairs(slopes_common, adjust = "tukey")

# Primary vs common-support treatment-specific slopes.
duration_sensitivity <- bind_rows(
  as.data.frame(slopes_primary) %>%
    mutate(Analysis = "Primary"),
  as.data.frame(slopes_common) %>%
    mutate(Analysis = paste0(
      "Common support (", common_min, "-", common_max, " days)"
    ))
)

print(duration_sensitivity)
write.csv(
  duration_sensitivity,
  file.path(tables_dir, "Sensitivity_common_duration_support.csv"),
  row.names = FALSE
)

# combo treatment
DiffplotCombo <- ggplot(
  data = WildLong2,
  aes(x = Dur, y = Diff, color = Combo, shape = Combo, fill = Combo)
) + 
  
  geom_smooth(
    method = "lm",
    se = TRUE,
    aes(color = Combo),
    alpha = 0.2
  ) + 
  
  # Redraw ONLY the control regression line on top
  geom_smooth(
    data = WildLong2[WildLong2$Combo == levels(WildLong2$Combo)[1], ],
    aes(x = Dur, y = Diff),
    method = "lm",
    se = FALSE,
    color = "gray40",
    linewidth = 1.2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  geom_point(
    aes(fill = Combo),
    size = 2,
    color = "black",
    show.legend = TRUE,
    position = position_jitter(width = 5)
  ) +
  
  theme_bw() + 
  
  labs(
    x = "Duration between Paired Samples", 
    y = expression(Delta * " Epigenetic Age (Days)")
  ) + 
  
  scale_shape_manual(
    values = c(21, 22, 24, 23),
    labels = c(
      "Control",
      "Anthelminthic",
      "Food",
      "Food + Anthelminthic"
    )
  ) + 
  
  scale_fill_manual(
    name = NULL, 
    values = c("gray50", "blue", "red", "purple"),
    labels = c(
      "Control",
      "Anthelminthic",
      "Food",
      "Food + Anthelminthic"
    )
  ) +
  
  scale_color_manual(
    name = NULL, 
    values = c("gray50", "blue", "red", "purple"),
    labels = c(
      "Control",
      "Anthelminthic",
      "Food",
      "Food + Anthelminthic"
    )
  ) +
  
  geom_abline(
    intercept = 0, 
    slope = 1, 
    size = 0.75, 
    color = "black", 
    linetype = "dashed"
  ) +
  
  guides(
    fill = guide_legend(title = NULL),
    color = guide_legend(title = NULL),
    shape = guide_legend(title = NULL)
  ) +
  
  theme(legend.position = "bottom")

DiffplotCombo

# Plot displayed in RStudio but not exported: EpiChangeGraphsCombo.pdf

DiffplotSex <- ggplot(data = WildLong2, aes(x = Sex, y = Diff, color = Repro.1)) + 
  geom_point(aes(shape = Sex), size = 2, show.legend = T) +
  theme_bw() + 
  labs(x = "Duration between Paired Samples", y = "Difference in EpiAge") + 
  scale_fill_manual(name = "", values = c("red","blue")) +
  scale_color_manual(name = "", values = c("red","blue"))
DiffplotSex

## option 2. do factors predict epigenetic ageing slopes?
longmod <- lm(Slope ~ Food + IV + Year + Sex + Repro.1, data = WildLong2)
summary(longmod)
#hist(resid(longmod))
anova(longmod)
nobs(longmod)

# slope graphs
slopegraph <- ggplot(WildLong2, aes(x=Food, y=Slope, fill=IV)) +
  geom_abline(intercept = 1, slope = 0, size = 1, color = "black", linetype = "dashed") + 
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = 0.15),aes(group=IV, size=Dur), col="black", pch=1) + 
  theme_bw() +
  labs(x="Food Supplementation", y = "Slope of EpiAge") + 
  scale_fill_manual(name = "", values = c("gray50","coral")) + 
  theme(text=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title = element_text(size=12),
        legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10))
slopegraph

reproslgraph <- ggplot(WildLong2, aes(x=Repro.1, y=Slope)) +
  geom_abline(intercept = 1, slope = 0, size = 1, color = "black", linetype = "dashed") + 
  geom_boxplot() +
  geom_point(col="black", pch=1, size = 2) + 
  theme_bw() +
  labs(x="Initial Reproductive Status", y = "Slope of Change in Epigenetic Age") + 
  theme(text=element_text(size=12), legend.text=element_text(size=12),
        axis.text=element_text(size=12), axis.title = element_text(size=12),
        legend.position = "bottom",
        legend.position.inside = c(0.3,0.1), plot.margin = ggplot2::margin(10,10,10,10))
reproslgraph

TmtSlopes <- ggplot(WildLong2, aes(x = Slope)) +   # map fill to Food here
  geom_histogram(bins = 50, color = "black", fill = "gray") +
  theme_bw() +
  facet_wrap(~Combo, nrow = 4) + 
  #scale_fill_manual(values = c("gray", "orange")) +   # custom colors for Food levels
  geom_vline(aes(xintercept = 1), color = "red", linewidth = 2, linetype = "dashed") + 
  theme(legend.position = c(0.8, 0.8),                    # fixed legend position
        legend.background = element_rect(color = "black"),
        legend.title = element_blank()) + 
  labs(x = "Slope of change in epigenetic age",
       y = "Number of Mice") + 
  theme(legend.position = "none")
TmtSlopes

AllSlopes <- ggplot(WildLong2, aes(x = Slope)) +   # map fill to Food here
  geom_histogram(bins = 50, color = "black", fill = "gray") + 
  theme_bw() + 
  geom_vline(aes(xintercept = 1), color = "black", linewidth = 2, linetype = "dashed") + 
  theme(legend.position = c(0.8, 0.8),                    # fixed legend position
        legend.background = element_rect(color = "black"),
        legend.title = element_blank()) + 
  labs(x = "Slope of Change in Epigenetic Age",
       y = "Number of Mice") + 
  theme(legend.position = "none")
AllSlopes

diff1 <- ggarrange(DiffplotIV, DiffplotFood, labels = c("B","C"), nrow = 1)


grDevices::pdf(file.path(figures_dir, "Figure_3.pdf"), width = 5, height = 6)
print(ggarrange(AllSlopes, DiffplotCombo, labels = c("A","B"), nrow = 2, heights = c(0.4,0.6)))
grDevices::dev.off()

# option 3. using pan mammalian epigenetic age as a proxy for biological age

PMmodM <- lm(predPMAge1 ~ predWMAge + Age2 + predWMAge*AnthelminthicTmt + predWMAge*SupplementTmt +
              Jul2 + BodyMass + Repro2 + Sex + Site + Year,
            data = Post, na.action = na.omit)
summary(PMmodM)
car::vif(PMmodM)
car::Anova(PMmodM, type = 2, test.statistic = "F")
nobs(PMmodM)
visreg(PMmodM, "SupplementTmt")
visreg(PMmodM, "AnthelminthicTmt")

Chron2 <- ggplot(
  data = Post,
  aes(x = scale(predWMAge), y = scale(predPMAge1), color = SupplementTmt, shape = SupplementTmt, fill = SupplementTmt)
) + 
  geom_point(size = 2, show.legend = TRUE) +
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() + 
  labs(
    x = "Wood Mouse Epigenetic Age", 
    y = "PM Clock #1 Epigenetic Age"
  ) + 
  scale_fill_manual(
    name = "Treatment", 
    values = c("gray", "blue", "red", "purple")
  ) +
  scale_color_manual(
    name = "Treatment", 
    values = c("gray", "blue", "red", "purple")
  ) +
  geom_abline(
    intercept = 0, 
    slope = 1, 
    size = 0.75, 
    color = "black", 
    linetype = "dashed"
  ) +
  guides(
    fill = guide_legend("Treatment"),
    color = guide_legend("Treatment"),
    shape = guide_legend("Treatment")
  )
Chron2

# Plot displayed in RStudio but not exported: Chron2GraphsCombo.pdf


## big epigAge figure
grDevices::pdf(file.path(figures_dir, "Figure_2.pdf"), width = 9, height = 6.2)
print(ggarrange(AgeCatJ, julgraph3,massgraph,ParNum,Ticks_plot,WMHV_plot,
          labels = c("A","B","C","D","E","F"), ncol = 3, nrow = 2))
grDevices::dev.off()


# Store the fully derived ecological data for use in downstream analyses.
# (This updates the intermediate snapshot saved following the initial data preparation.)
save(data, file = file.path(processed_dir, "WMClock_EcologicalData.Rdata"))

# Save longitudinal summaries and sensitivity subsets together for reuse.
save(WildLong, WildLong2, WildLong_no_extreme, WildLong_common,
     file = file.path(processed_dir, "WMClock_LongitudinalData.Rdata"))

# End of Script 04. Only manuscript figures are exported.
message("Script 04 complete. Outputs: ", results_dir)
