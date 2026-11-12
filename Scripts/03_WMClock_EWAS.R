################################################################################
# 03_WMClock_EWAS.R
# Age-associated DNA methylation in laboratory wood mice
################################################################################
# PURPOSE
#   Identify CpG sites associated with chronological age in laboratory wood
#   mouse ear samples, adjusting for sex. Correct EWAS test-statistic inflation
#   using bacon, annotate age-associated CpGs, and characterise enriched
#   functional categories using GREAT (human hg19 genomic coordinates).
#
# RUN FROM PROJECT ROOT: WoodMouseEpigeneticClock/
# INPUT FILES
#   Data/Processed Data/beta_pheno_anno.Rdata
#       RData object joined2mcpg: batch-corrected beta values and sample data
#   Data/PanMamManifest_Aligned_ASylvaticus.csv
#       Wood mouse CpG annotations (CGid, geneChr, geneStart, geneEnd, etc.)
#   Data/Manifest, HorvathMammalMethylChip40.csv
#       Array manifest with human hg19 mapping (online GREAT mode only)
#   Results/03_WMClock_EWAS/Tables/GREATResults.csv
#   Results/03_WMClock_EWAS/Tables/GREATResults_Top.csv
#       Archived enrichment tables (required when run_great_online = FALSE)
#   Data/Processed Data/WMClock_Figure1_plots.rds
#       Figure 1A-C ggplot objects saved by Script 02
#
# OUTPUTS
#   Results/03_WMClock_EWAS/Tables/
#       EWAS_sample_summary_by_batch.csv
#       EWAS_genomic_inflation.csv
#       EWAS_age_bacon_adjusted.csv
#       GREATResults.csv
#       GREATResults_Top.csv
#   Results/03_WMClock_EWAS/Figures/
#       EWAS_age_QQ_plots.pdf
#       Figure_1_panels_D-E.pdf
#       Figure_1_panels_F-G.pdf
#       Figure_1_panels_D-G.pdf
#       Figure_1_Complete.pdf  (combined publication figure A-G)
#
# ANALYSIS NOTES
#   * Age effect: limma linear model ~ Age.Days + Sex, fitted per CpG.
#   * bacon adjusts EWAS t-statistics and p-values; BH FDR controls multiplicity.
#   * GREAT uses human hg19 orthologous mapping, analysed separately for CpGs
#     positively and negatively associated with age.
#   * Figure 1A-C: clock training/test and pan-mammalian validation plots
#     loaded from the saved ggplot objects created by Script 02.
#   * Figure 1D/E: Manhattan plot and genomic annotation categories.
#     Figure 1F/G: GREAT enrichment for hyper/hypomethylated CpGs.
#   * Set run_great_online = FALSE to use archived GREAT tables while the server
#     is unavailable. Set TRUE to rerun GREAT using the original settings.
#   * Archived enrichment results reflect the EWAS run that originally generated
#     those files; they are not recomputed from the current EWAS results.
################################################################################

# ---- Analysis switch ----------------------------------------------------------
# FALSE: skip the GREAT server; use archived GREATResults.csv and
#        GREATResults_Top.csv from this analysis's Tables folder.
# TRUE:  rerun GREAT on the significant CpGs and overwrite those tables.
run_great_online <- FALSE

# ---- 1. Packages, project paths, and input validation ------------------------
suppressPackageStartupMessages({
  library(limma)
  library(qqman)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(bacon)
  # rGREAT is loaded only when online GREAT submission is enabled.
})

# Use relative paths with the R project root as the working directory.
input_file <- file.path("Data", "Processed Data", "beta_pheno_anno.Rdata")
clock_plots_file <- file.path("Data", "Processed Data", "WMClock_Figure1_plots.rds")
annotation_file <- file.path("Data", "PanMamManifest_Aligned_ASylvaticus.csv")
manifest_file <- file.path("Data", "Manifest, HorvathMammalMethylChip40.csv")
result_dir <- file.path("Results", "03_WMClock_EWAS")
table_dir <- file.path(result_dir, "Tables")
figure_dir <- file.path(result_dir, "Figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(input_file, annotation_file, clock_plots_file)
if (run_great_online) {
  required_files <- c(required_files, manifest_file)
} else {
  required_files <- c(required_files,
                      file.path(table_dir, "GREATResults.csv"),
                      file.path(table_dir, "GREATResults_Top.csv"))
}
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

# ---- 2. Select known-aged laboratory ear samples -----------------------------
# ---- load methylation data ---------------------------------------------------
load(input_file)
stopifnot(exists("joined2mcpg"), is.data.frame(joined2mcpg))
joined2mcpg <- joined2mcpg[!is.na(joined2mcpg$Age.Days) & joined2mcpg$Where == "Lab" & joined2mcpg$Tissue == "Ear",]

summary_by_batch <- joined2mcpg %>%
  group_by(Batch) %>%
  dplyr::summarise(
    n_rows    = n(),
    age_avg   = mean(Age.Days, na.rm = TRUE),
    age_min   = min(Age.Days, na.rm = TRUE),
    age_max   = max(Age.Days, na.rm = TRUE),
    age_range = max(Age.Days, na.rm = TRUE) - min(Age.Days, na.rm = TRUE),
    age_se    = sd(Age.Days, na.rm = TRUE) / sqrt(sum(!is.na(Age.Days))),
    n_male    = sum(tolower(as.character(Sex)) %in% c("m", "male"), na.rm = TRUE),
    n_female  = sum(tolower(as.character(Sex)) %in% c("f", "female"), na.rm = TRUE),
    .groups = "drop"
  )

summary_by_batch
write.csv(summary_by_batch, file.path(table_dir, "EWAS_sample_summary_by_batch.csv"), row.names = FALSE)

# ---- 3. Prepare CpG-by-sample matrix and limma design -------------------------

cg_cols <- grep("cg", colnames(joined2mcpg), ignore.case = TRUE)
# Sample identifier is selected by name to avoid assuming a fixed metadata order.
betas <- joined2mcpg[ , c("Basename", colnames(joined2mcpg)[cg_cols])]
rownames(betas) <- betas$Basename; betas <- as.data.frame(t(betas[,-1]))

pheno <- joined2mcpg[ ,(!grepl("cg", colnames(joined2mcpg), ignore.case = TRUE))]

# ---- prepare data ------------------------------------------------------------
# betas: CpGs × samples, already batch-corrected.
# pheno: matching sample metadata, in the same order as beta columns.

# Check order
stopifnot(all(pheno$Basename == colnames(betas)))

# Build design matrix (test association with age, control for sex)
design <- model.matrix(~ Age.Days + Sex, data = pheno)

# ---- 4. Fit age EWAS and correct inflation with bacon --------------------------
fit <- lmFit(betas, design)
fit2 <- eBayes(fit)

# Extract t-statistics and p-values for Age
tvals_age <- fit2$t[ , "Age.Days"]
pvals_age <- fit2$p.value[ , "Age.Days"]
cor <- fit2$coefficients[ , "Age.Days"]

# ---- apply bacon to correct inflation ----------------------------------------
# bacon expects t-statistics as input
bc <- bacon(teststatistics = tvals_age)

# Extract bacon-adjusted t-statistics and p-values
tvals_bacon <- bacon::tstat(bc)
pvals_bacon <- pval(bc)

# ---- Measure genomic inflation before and after bacon -------------------------
# before bacon
chisq_before <- qchisq(1 - pvals_age, df=1)
lambda_before <- median(chisq_before, na.rm=TRUE) / qchisq(0.5, df=1)

# after bacon
chisq_after <- qchisq(1 - pvals_bacon, df=1)
lambda_after <- median(chisq_after, na.rm=TRUE) / qchisq(0.5, df=1)

cat("Genomic inflation factor (λ) before bacon: ", round(lambda_before,3), "\n")
cat("Genomic inflation factor (λ) after  bacon: ", round(lambda_after,3), "\n")
write.csv(data.frame(Stage = c("Before bacon", "After bacon"),
                     Lambda = c(lambda_before, lambda_after)),
          file.path(table_dir, "EWAS_genomic_inflation.csv"), row.names = FALSE)

# ---- QQ plots ----------------------------------------------------------------
# Save the two base-R QQ plots together without modifying their style.
pdf(file.path(figure_dir, "EWAS_age_QQ_plots.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
# QQ plot before bacon
qq(pvals_age, main = "QQ plot: EWAS for age (before bacon)")

# QQ plot after bacon
qq(pvals_bacon, main = "QQ plot: EWAS for age (after bacon)")
dev.off()

# ---- 5. Annotate age-associated CpGs and save EWAS table ----------------------
# save adjusted p-values and t-values
EWASResults <- data.frame(
  CpG = rownames(betas),
  t_bacon = tvals_bacon,
  p_bacon = pvals_bacon,
  cor
)

# fdr-adjusted pvalues
EWASResults$pvaladj <- p.adjust(EWASResults$p_bacon, method = "fdr")

# number of significant CpGs associated with age
nrow(EWASResults[EWASResults$pvaladj<0.05 & EWASResults$cor < 0,])


# ---- load and prepare annotated CpG information ------------------------------
cpgs <- read.csv(annotation_file, stringsAsFactors = FALSE)

cpgsM <- cpgs %>%
  dplyr::select(
    CPG = CGid,
    MChr = geneChr,
    MStart = geneStart,
    MEnd = geneEnd,
    MGene = geneId,
    MCat = annotation2
  ) %>%
  mutate(
    MCat = factor(MCat, levels = c("Downstream", "3' UTR", "Intron", "Exon", "5' UTR", "Promoter", "Intergenic")),
    MChr = as.factor(MChr)
  )

# ---- merge EWAS results with annotation --------------------------------------
results_bacon <- EWASResults %>%
  left_join(cpgsM, by = c("CpG" = "CPG")) 

results_bacon$MChr <- as.factor(results_bacon$MChr)
results_bacon$MChr <- factor(results_bacon$MChr, levels=c("1"="1","2"="2","3"="3","4"="4","5"="5","6"="6",
                                                          "7"="7","8"="8","9"="9","10"="10","11"="11",
                                                          "12"="12","13"="13","14"="14","15"="15",
                                                          "16"="16","17"="17","18"="18","19"="19",
                                                          "20"="20","21"="21","22"="22","23"="23","24"="X","29"="Y"))
results_bacon$ChrCol <- as.character(plyr::revalue(results_bacon$MChr, c("1" = "gray","2" = "black","3" = "gray","4" = "black","5" = "gray",
                                                                   "6" = "black","7" = "gray","8" = "black","9" = "gray","10" = "black",
                                                                   "11" = "gray","12" = "black","13" = "gray","14" = "black","15" = "gray",
                                                                   "16" = "black","17" = "gray","18" = "black","19" = "gray","20" = "black",
                                                                   "21" = "gray","22" = "black","23" = "gray","24"="gray","X" = "black","Y" = "gray")))
results_bacon$Direction <- as.factor(ifelse(results_bacon$cor > 0, "Positive", "Negative"))
results_bacon$ChrCol[which(results_bacon$Direction == "Positive" & results_bacon$pvaladj < 0.05)] <- "red"
results_bacon$ChrCol[which(results_bacon$Direction == "Negative" & results_bacon$pvaladj < 0.05)] <- "blue"
results_bacon$ChrCol <- as.factor(results_bacon$ChrCol)
results_bacon$MGene <- toupper(results_bacon$MGene)
results_bacon <- results_bacon[order(results_bacon$MChr, results_bacon$MStart), ]
results_bacon$order <- seq.int(nrow(results_bacon))
results_bacon$log10P <- -1 * log10(results_bacon$pvaladj)

write.csv(results_bacon, file.path(table_dir, "EWAS_age_bacon_adjusted.csv"), row.names=FALSE)

# ---- 6. Make Figure 1D/E: Manhattan plot and CpG categories -------------------
top_genes <- results_bacon %>%
  group_by(MGene) %>%
  slice_max(log10P, n = 1, with_ties = FALSE) %>%  # top site per gene
  ungroup() %>%
  slice_max(log10P, n = 50)  # top 50 overall

manhattan_plot <- ggplot(results_bacon, aes(x = order, y = log10P, color = ChrCol, label = MGene)) +
  geom_point(size = 2) +
  labs(x = "Chromosome", y = "-log10(P)") +
  theme_classic() +
  scale_color_manual(values = c("black", "blue", "gray", "red")) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.08, "cm")
  ) +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  geom_label_repel(
    data = top_genes,
    label.size = NA, fill = alpha("white", 0), color = "black",
    size = 3, max.overlaps = 15
  )
manhattan_plot

# ---- graph the CpG distribution ----------------------------------------------
# barplot: significant CpG locations
EWAS_cpg_cat <- results_bacon %>%
  filter(pvaladj <= 0.05) %>%
  ggplot(aes(x = MCat, fill = Direction)) +
  geom_bar(position = "dodge", color = "black") +
  coord_flip() +
  xlab("Closest TSS to CpG") +
  ylab("") +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    text = element_text(size = 10)
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(values = c("blue", "red"), labels = c("hypomethylated", "hypermethylated")) +
  theme(legend.position = "bottom", legend.direction = "horizontal")

# barplot: background CpG distribution
All_cpg_cat_background <- ggplot(cpgsM, aes(x = MCat)) +
  geom_bar(color = "black") +
  coord_flip() +
  xlab("CpG Location Relative to Closest TSS") +
  ylab("") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.5, 0.1),
    legend.title = element_blank(),
    text = element_text(size = 10)
  ) +
  scale_x_discrete(drop = FALSE)

# arrange both plots side by side
rec1 <- ggarrange(EWAS_cpg_cat, All_cpg_cat_background, nrow = 1,
                  widths = c(0.6, 0.4), hjust = c(-0.5, 0.7), common.legend = TRUE)

# Preserve the original plotting layers, arrangement, panel widths and size.
# Manuscript panel letters are D and E.
thirdrow <- ggarrange(manhattan_plot, rec1,
                      labels = c("D", "E"), ncol = 2, widths = c(0.6, 0.4))
pdf(file.path(figure_dir, "Figure_1_panels_D-E.pdf"), width = 8.19, height = 3.92)
print(thirdrow)
dev.off()

# ---- 7. GREAT functional enrichment (human hg19 coordinates) -----------------
# Significant CpGs are joined to the mammalian array manifest. The input must
# provide the same human coordinate fields used for the GREAT intervals.
# GREAT is run separately for sites whose methylation increases vs decreases
# with chronological age.

if (run_great_online) {
  # Reproduce the GREAT 3 analysis using the human hg19 mappings. These
  # requests require a functioning GREAT server and will overwrite archived
  # enrichment tables in this folder when successful.
  library(rGREAT)
  mani <- read.csv(manifest_file, header = TRUE, stringsAsFactors = TRUE)

  # ---- isolate datasets of significant CpGs that either increase or decrease with age
  diffH <- left_join(results_bacon,mani,by=c("CpG" = "IlmnID"))
  diffH <- diffH %>% subset(pvaladj < 0.05) 

  # Increase with age: prepare genomic intervals for GREAT.
  # These numeric column positions are based on the input manifest structure;
  # check their identities if the manifest format changes.
  diffI <- diffH[diffH$Direction == "Positive",]
  diffI <- diffI[,c(44,41,42,1)]
  diffI <- diffI[diffI$Human.Hg19_Chr!="",]
  colnames(diffI) <- c("chr","start","end","cpg")
  diffI$chr <- droplevels(diffI$chr)

  # decrease with age
  diffD <- diffH[diffH$Direction=="Negative",]
  diffD <- diffD[,c(44,41,42,1)]
  diffD <- diffD[diffD$Human.Hg19_Chr!="",]
  colnames(diffD) <- c("chr","start","end","cpg")
  diffD$chr <- droplevels(diffD$chr)

  # perform the GREAT analysis for significantly increasing CpGs
  output.all={}
  job1 = submitGreatJob(diffI,
                       species               = "hg19",
                       includeCuratedRegDoms = TRUE,
                       rule                  = c("basalPlusExt"),
                       adv_upstream          = 5.0,
                       adv_downstream        = 1.0,
                       adv_span              = 50,
                       request_interval = 300,
                       version="3",
                       max_tries = 10)

  ontology.all=availableOntologies(job1)
  print(ontology.all)

  for(k in 1:length(ontology.all)){
    print(ontology.all[k])
    out0.list = tryCatch(getEnrichmentTables(job1,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
    if(!is.null(out0.list)){
      db0.list=as.list(names(out0.list))
      output<-Map(cbind,Database=db0.list,out0.list)
      output<-do.call('rbind',output)
      output.all=rbind(output.all,output)
    }
  }

  # save each ontology output separately and then rbind together
  MouseGMF <- getEnrichmentTable(job1, ontology = "GO Molecular Function", download_by = "tsv")
  MouseGBP <- getEnrichmentTable(job1, ontology = "GO Biological Process", download_by = "tsv")
  MouseGCC <- getEnrichmentTable(job1, ontology = "GO Cellular Component", download_by = "tsv")
  MouseMF <- getEnrichmentTable(job1, ontology = "Mouse Phenotype", download_by = "tsv")
  MouseHF <- getEnrichmentTable(job1, ontology = "Human Phenotype", download_by = "tsv")
  MouseDO <- getEnrichmentTable(job1, ontology = "Disease Ontology", download_by = "tsv")
  MousePP <- getEnrichmentTable(job1, ontology = "PANTHER Pathway", download_by = "tsv")
  MouseSG <- getEnrichmentTable(job1, ontology = "MSigDB Perturbation", download_by = "tsv")

  names(MouseSG)[2] <- "Desc"
  names(MouseSG)[3] <- "ID"

  MouseAll1 <- rbind(MouseGMF,MouseGBP,MouseGCC,MouseMF,MouseHF,MouseDO,MousePP,MouseSG)
  MouseAll1$Direction <- "Increase"

  # perform the GREAT analysis for significantly decreasing CpGs
  output.all={}
  job2 = submitGreatJob(diffD,
                       species               = "hg19",
                       includeCuratedRegDoms = TRUE,
                       rule                  = c("basalPlusExt"),
                       adv_upstream          = 5.0,
                       adv_downstream        = 1.0,
                       adv_span              = 50,
                       request_interval = 300,
                       version="3",
                       max_tries = 10)

  ontology.all=availableOntologies(job2)
  print(ontology.all)

  for(k in 1:length(ontology.all)){
    print(ontology.all[k])
    out0.list = tryCatch(getEnrichmentTables(job2,ontology=ontology.all[k], download_by = "tsv"),error=function(e){NULL})
    if(!is.null(out0.list)){
      db0.list=as.list(names(out0.list))
      output<-Map(cbind,Database=db0.list,out0.list)
      output<-do.call('rbind',output)
      output.all=rbind(output.all,output)
    }
  }

  # save each ontology output separately and then rbind together
  MouseGMF <- getEnrichmentTable(job2, ontology = "GO Molecular Function", download_by = "tsv")
  MouseGBP <- getEnrichmentTable(job2, ontology = "GO Biological Process", download_by = "tsv")
  MouseGCC <- getEnrichmentTable(job2, ontology = "GO Cellular Component", download_by = "tsv")
  MouseMF <- getEnrichmentTable(job2, ontology = "Mouse Phenotype", download_by = "tsv")
  MouseHF <- getEnrichmentTable(job2, ontology = "Human Phenotype", download_by = "tsv")
  MouseDO <- getEnrichmentTable(job2, ontology = "Disease Ontology", download_by = "tsv")
  MousePP <- getEnrichmentTable(job2, ontology = "PANTHER Pathway", download_by = "tsv")
  MouseSG <- getEnrichmentTable(job2, ontology = "MSigDB Perturbation", download_by = "tsv")

  names(MouseSG)[2] <- "Desc"
  names(MouseSG)[3] <- "ID"

  MouseAll2 <- rbind(MouseGMF,MouseGBP,MouseGCC,MouseMF,MouseHF,MouseDO,MousePP,MouseSG)
  MouseAll2$Direction <- "Decrease"


  # add direction label
  MouseAll1$Direction <- "Hypermethylated (↑ with age)"
  MouseAll2$Direction <- "Hypomethylated (↓ with age)"

  # combine both datasets
  great_combined <- bind_rows(MouseAll1, MouseAll2)
  write.csv(great_combined, file.path(table_dir, "GREATResults.csv"), row.names = FALSE)

  GREAT <- read.csv(file.path(table_dir, "GREATResults.csv"), header = TRUE, stringsAsFactors = TRUE)

  # filter for significant terms
  great_combined <- GREAT %>%
    filter(BinomFdrQ < 0.05) %>%
    mutate(
      log10P = -log10(BinomFdrQ),
      Desc = fct_reorder(Desc, log10P)  # order within plot
    )

  # for each direction and ontology set, select top 10 terms
  great_top <- great_combined %>%
    group_by(Direction, Ontology) %>%
    slice_max(order_by = log10P, n = 10) %>%
    ungroup() %>%
    mutate(
      Term = str_trunc(as.character(Desc), 30),   
      Term = fct_reorder(Term, GeneFoldEnrich) 
    )

  great_top$Term <- str_wrap(great_top$Term, width = 50)

  write.csv(great_top, file.path(table_dir, "GREATResults_Top.csv"), row.names = FALSE)
  great_top <- read.csv(file.path(table_dir, "GREATResults_Top.csv"), header = TRUE, stringsAsFactors = TRUE)


} else {
  # GREAT web service is unavailable: use the archived enrichment outputs.
  # Both are pre-existing results, not newly calculated by this execution.
  GREAT <- read.csv(file.path(table_dir, "GREATResults.csv"),
                    header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE)
  great_top <- read.csv(file.path(table_dir, "GREATResults_Top.csv"),
                        header = TRUE, stringsAsFactors = FALSE,
                        check.names = FALSE)

  # Confirm that the saved files contain the fields needed downstream.
  stopifnot(all(c("Ontology", "Desc", "BinomFdrQ", "Direction") %in% names(GREAT)))
  stopifnot(all(c("Ontology", "Desc", "Direction", "log10P") %in% names(great_top)))
  message("Using archived GREAT enrichment results from: ", table_dir)
  message("GREAT web submissions skipped (run_great_online = FALSE).")
}

# ---- 8. Make Figure 1F/G: GREAT enrichment plots ------------------------------
# Retain the original two ontology panels, scales, labels and plotting style.
# graph of enrichment for hypermethylated CpGs
pI <- ggplot(great_top[great_top$Direction=='Hypermethylated (↑ with age)' & great_top$Ontology %in% c('GO Molecular Function','GO Biological Process'),], aes(x = log10P, y = Desc)) +
  geom_segment(aes(x = 0, y = Desc, xend = log10P, yend = Desc), color = "gray60") +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  coord_cartesian(clip = "off") +
  facet_wrap(~ Ontology, scales = "free_y", nrow = 2) +
  labs(
    x = "log10(FDR p-value)",
    y = NULL
  ) +
  xlim(0,250) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    #plot.margin = margin(t = 5, r = 5, b = 5, l = 55),
    axis.text.y = element_text(size = 9)
  )
pI


# graph of enrichment for hypomethylated CpGs
pD <- ggplot(great_top[great_top$Direction=='Hypomethylated (↓ with age)' & great_top$Ontology %in% c('GO Molecular Function','GO Biological Process'),], aes(x = log10P, y = Desc)) +
  geom_segment(aes(x = 0, y = Desc, xend = log10P, yend = Desc), color = "gray60") +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  coord_cartesian(clip = "off") +
  facet_wrap(~ Ontology, scales = "free_y", nrow = 2) +
  labs(
    x = "log10(FDR p-value)",
    y = NULL
  ) +
  xlim(0,250) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    #plot.margin = margin(t = 5, r = 5, b = 5, l = 100),
    axis.text.y = element_text(size = 9)
  )
pD

# Save the GREAT panels with the original panel labels and dimensions.
secrow <- ggarrange(pI, pD,
                    labels = c("F. Hypermethylated", "G. Hypomethylated"), ncol = 2)
pdf(file.path(figure_dir, "Figure_1_panels_F-G.pdf"), width = 12, height = 5)
print(secrow)
dev.off()

# Standalone collection of this script's panels D-G. Figure 1A-C is created
# by the clock-creation workflow and is not required to run this EWAS script.
pdf(file.path(figure_dir, "Figure_1_panels_D-G.pdf"), width = 11, height = 8)
print(ggarrange(thirdrow, secrow, ncol = 1))
dev.off()


# ---- 9. Assemble complete manuscript Figure 1 (panels A-G) -------------------
# Import native ggplot objects from the clock-creation workflow. Their original
# point geometry, themes, fitted lines, axes, and automated accuracy labels are
# preserved; no PDF importing, rasterisation, or clock-model refitting is needed.
clock_plots <- readRDS(clock_plots_file)
needed_plot_names <- c("All_graph", "TestAll", "PanOnly_Graph2")
if (!is.list(clock_plots) ||
    !all(needed_plot_names %in% names(clock_plots)) ||
    !all(vapply(clock_plots[needed_plot_names],
                function(p) inherits(p, "ggplot"), logical(1)))) {
  stop("The saved clock plot file must contain ggplot objects named: ",
       paste(needed_plot_names, collapse = ", "),
       ". Rerun Script 02 to regenerate it.")
}

# Top row: original clock-validation plots (A-C).
firstrow <- ggarrange(
  clock_plots$All_graph,
  clock_plots$TestAll,
  clock_plots$PanOnly_Graph2,
  labels = c("A", "B", "C"),
  ncol = 3
)

# The EWAS row (D-E) is stored in 'thirdrow' above and the GREAT row (F-G)
# in 'secrow'. Reuse them without changing their plot aesthetics or legends.
Figure1 <- ggarrange(
  firstrow,
  thirdrow,
  secrow,
  ncol = 1,
  heights = c(1, 1, 1)
)

# Opening a PDF graphics device and printing ggplot objects keeps text, lines,
# points, and axis elements vector-based, suitable for journal submission.
figure1_output <- file.path(figure_dir, "Figure_1_Complete.pdf")
pdf(figure1_output, width = 11, height = 11)
print(Figure1)
dev.off()
message("Saved complete vector Figure 1 (A-G): ", figure1_output)

message("EWAS and GREAT figures complete (archived GREAT tables used when online mode is off).\nTables: ", table_dir,
        "\nFigures: ", figure_dir)
