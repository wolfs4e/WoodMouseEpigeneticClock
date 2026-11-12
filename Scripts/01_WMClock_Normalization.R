################################################################################
# 01_WMClock_Normalization.R
# Wood mouse epigenetic clock | IDAT QC, normalization, batch correction, joining
################################################################################
# PURPOSE
#   Process MammalMethylChip40 IDAT files from three datasets (2023, original
#   colony dataset, and 2024). Apply detection-p quality control, normalize
#   using minfi's Noob method, remove probes failing QC in any dataset, and
#   correct for batch effects using ComBat. Finally, join the corrected beta
#   values to sample metadata and retain CpGs represented in the wood mouse
#   annotation/manifest file.
#
# PROJECT STRUCTURE (run from the R project root / GitHub repository root)
#   Data/
#     DNAm 2023/Idats/
#     DNAm 2024/Idats2/
#     DNAm Original Colony Dataset/idatFiles/
#     WMClock_3Batches_SampleInfo.csv
#     PanMamManifest_Aligned_ASylvaticus.csv
#     Processed Data/
#   Scripts/01_WMClock_Normalization.R
#   Results/
#
# INPUTS
#   1. Three raw IDAT directories.
#   2. Data/WMClock_3Batches_SampleInfo.csv: sample phenotype/batch metadata.
#   3. Data/PanMamManifest_Aligned_ASylvaticus.csv: CpG annotation manifest.
#
# FINAL OUTPUTS
#   Data/Processed Data/beta_pheno_anno.Rdata
#   Data/Processed Data/beta_pheno_anno.csv
#   Both contain joined2mcpg: sample metadata + retained annotated CpG betas.
#
################################################################################

# ---- 1. Load packages --------------------------------------------------------
suppressPackageStartupMessages({
  library(minfi)
  library(magrittr)
  library(bacon)
  library(sva)
  library(dplyr)
})

# ---- 2. Input paths and quality-control thresholds ---------------------------
# The three input datasets are normalized separately before merging.
base_dirs <- c("Data/DNAm 2023/Idats",
               "Data/DNAm Original Colony Dataset/idatFiles",
               "Data/DNAm 2024/Idats2")

# A sample-probe measurement fails when detection P > 0.01.
# Remove samples when >20% of probes fail; flag probes when >10% of remaining
# samples fail. Flagged CpGs are removed from ALL datasets after processing.
p_thresh        <- 0.01   # CpG-sample detection-P failure cut-off
max_fail_sample <- 0.20   # maximum fraction of failed CpGs per sample
max_fail_probe  <- 0.10   # maximum fraction of failed samples per CpG

# ---- 3. Read, quality-control and Noob-normalize each IDAT dataset ------------
process_dataset <- function(path) {
  message("\nProcessing: ", path)

  # Read the sample sheet and corresponding raw red/green intensities.
  targets <- read.metharray.sheet(path, recursive = TRUE)
  targets$Basename2 <- paste0(targets$Slide, "_", targets$Array)

  RG <- read.metharray.exp(base = NULL, targets = targets, recursive = TRUE)

  # Identify the mammalian methylation array to minfi.
  RG@annotation <- c(array = "HorvathMammalMethylChip40",
                     annotation = "test.unknown")

  # Calculate detection P for every CpG-sample measurement.
  detP  <- detectionP(RG)
  fail  <- detP > p_thresh

  # Exclude samples with more than 20% failed CpGs.
  bad_samp <- colMeans(fail) > max_fail_sample
  message("  • removing ", sum(bad_samp), " bad samples (",
          round(100 * mean(bad_samp), 1), "%)")

  RG    <- RG[ , !bad_samp]
  detP  <- detP[ , !bad_samp]
  fail  <- fail  [ , !bad_samp]

  # Flag CpGs failing in >10% of the retained samples in THIS dataset.
  # Do not drop them here; all datasets will share one combined exclusion list.
  bad_cpg <- rowMeans(fail) > max_fail_probe
  message("  • flagging ", sum(bad_cpg), " bad CpGs (",
          round(100 * mean(bad_cpg), 1), "%)")

  # Noob background correction and dye-bias normalization; extract beta values.
  # The beta matrix is oriented CpGs (rows) x samples (columns).
  M.noob <- preprocessNoob(RG)
  beta   <- getBeta(M.noob)

  # Return normalized betas plus IDs/count of flagged CpGs for later filtering.
  list(beta      = beta,
       bad_cpg   = names(which(bad_cpg)),
       n_bad_cpg = sum(bad_cpg))
}

# Run the same preprocessing function independently on all three datasets.
results <- lapply(base_dirs, process_dataset)
names(results) <- basename(base_dirs)

# Print the number of CpGs flagged within each dataset.
for (i in seq_along(results)) {
  message("Dataset ", names(results)[i], ": flagged ",
          results[[i]]$n_bad_cpg, " CpGs for removal.")
}

# ---- 4. Apply consistent probe filtering and merge datasets ------------------
# Union: a CpG flagged in ANY dataset is removed from EVERY dataset.
union_bad_cpg <- unique(unlist(lapply(results, `[[`, "bad_cpg")))
message("\nTotal CpGs to remove across ALL datasets: ",
        length(union_bad_cpg))

# Remove the same union of flagged CpGs from each normalized beta matrix.
beta_final <- lapply(results, function(res) {
  res$beta[!rownames(res$beta) %in% union_bad_cpg, ]
})

# Keep CpGs shared by all three matrices, then concatenate sample columns.
common_cpgs <- Reduce(intersect, lapply(beta_final, rownames))
beta_merged <- do.call(cbind, lapply(beta_final, `[`, common_cpgs, ))

# ---- 5. Match sample metadata and prepare ComBat covariates -------------------
# Sample sheet includes Basename, Sex, Age.Days and Batch.
pheno <- read.csv("Data/WMClock_3Batches_SampleInfo.csv", header=T, stringsAsFactors = T)

# Ensure sample identifiers are character strings in both datasets.
colnames(beta_merged) <- as.character(colnames(beta_merged))
pheno$Basename <- as.character(pheno$Basename)

# Only retain metadata for samples that survived IDAT-level QC.
pheno_sub <- pheno %>% filter(Basename %in% colnames(beta_merged))

# Samples lacking a recorded sex are excluded from batch correction and output.
pheno_sub$Sex <- as.character(pheno_sub$Sex)
pheno_sub_clean <- pheno_sub %>%
  filter(!is.na(Sex) & Sex != "")

# Reorder methylation columns to the phenotype sample order; verify alignment.
beta_merged_clean <- beta_merged[, pheno_sub_clean$Basename]

pheno_sub_clean <- pheno_sub_clean[match(colnames(beta_merged_clean), pheno_sub_clean$Basename), ]
stopifnot(all(pheno_sub_clean$Basename == colnames(beta_merged_clean)))

# Missing ages are filled with the mean of observed ages for the ComBat design.
# Age.Days itself is not overwritten; the imputed values are in Age_imputed.
pheno_sub_clean$Age_imputed <- pheno_sub_clean$Age.Days
mean_age <- mean(pheno_sub_clean$Age.Days, na.rm = TRUE)
pheno_sub_clean$Age_imputed[is.na(pheno_sub_clean$Age_imputed)] <- mean_age

# Treat sex and processing batch as categorical variables.
pheno_sub_clean$Sex <- factor(pheno_sub_clean$Sex)
pheno_sub_clean$Batch <- factor(pheno_sub_clean$Batch)

# Include age and sex in the ComBat model so their associated variation is
# preserved during adjustment for the technical batch variable.
mod <- model.matrix(~ Age_imputed + Sex, data = pheno_sub_clean)

# ---- 6. Correct technical batch effects with ComBat --------------------------
# ComBat is applied to logit-transformed methylation beta values (M values).
mvals <- log2(beta_merged_clean / (1 - beta_merged_clean))

mvals_combat <- ComBat(
  dat = as.matrix(mvals),
  batch = pheno_sub_clean$Batch,
  mod = mod,
  par.prior = TRUE,
  prior.plot = FALSE
)

# Transform batch-corrected M values back to methylation proportions (betas).
beta_combat <- 2^mvals_combat / (1 + 2^mvals_combat)
message("Batch correction complete. Samples with missing sex were removed.")

# ---- 7. Join corrected methylation values to sample metadata -----------------
# Transpose to samples (rows) x CpGs (columns) for one row per sample.
dnam <- as.data.frame(t(beta_combat))
dnam$joinkey <- rownames(dnam)

# Join methylation columns to their matching sample/phenotype information.
joined <- left_join(pheno_sub_clean, dnam, by = c("Basename" = "joinkey"))

# ---- 8. Identify CpGs covered by the wood mouse manifest ---------------------
# Load the A. sylvaticus CpG manifest and construct an annotation lookup.
cpgs <- read.csv("Data/PanMamManifest_Aligned_ASylvaticus.csv", header=TRUE, stringsAsFactors = T)
cpgs$Chr <- as.factor(cpgs$geneChr)
cpgsM <- cpgs[,c(2,13,14,15,18,12)]
colnames(cpgsM) <- c("CPG","MChr","MStart","MEnd","MGene","MCat")
cpgsM$MCat <- factor(cpgsM$MCat, levels = c("Downstream","3' UTR","Intron","Exon","5' UTR","Promoter","Intergenic"))

# Preserve methylation beta columns whose CpG IDs are present in the manifest.
joined2mcpg <- cbind(joined[, 2:14],
  joined[, colnames(joined) %in% cpgsM[,1]])

# ---- 9. Save ONLY the final joined and filtered dataset -----------------------
dir.create("Data/Processed Data", recursive = TRUE, showWarnings = FALSE)
save(joined2mcpg, file = "Data/Processed Data/beta_pheno_anno.Rdata")
#write.csv(joined2mcpg, file = "Data/Processed Data/beta_pheno_anno.csv")

message("Final dataset saved in Data/Processed Data/ (Rdata and CSV).")
