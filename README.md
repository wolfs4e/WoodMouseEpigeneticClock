# Wood Mouse Epigenetic Clock

### DNA methylation ageing in laboratory and wild wood mice (*Apodemus sylvaticus*)

This repository contains the R scripts, processed datasets, and analytical outputs associated with the development and validation of a species-specific epigenetic clock for the wood mouse (*Apodemus sylvaticus*), and its application to investigate epigenetic ageing in wild populations.

The project combines DNA methylation data from laboratory and wild wood mice with longitudinal ecological observations to investigate age-associated methylation patterns, individual variation in epigenetic ageing, and the potential effects of environmental conditions on ageing trajectories.

## Research objectives

The analyses address four primary questions:

1. **Epigenetic clock development:** Can DNA methylation accurately predict chronological age in wood mice?
2. **Age-associated methylation:** Which CpG sites and biological pathways exhibit age-associated changes in DNA methylation?
3. **Ecological correlates of epigenetic age:** How does epigenetic age vary with sex, body condition, reproduction, seasonality, and parasite infection in wild populations?
4. **Environmental effects on ageing:** Do experimental food supplementation and anthelminthic treatment influence longitudinal epigenetic ageing trajectories in wild wood mice?

## Study system

The project uses wood mice sampled from laboratory colonies and wild populations in Scotland.

DNA methylation was measured using the Horvath Mammalian Methylation Array (MammalMethylChip40). A species-specific epigenetic clock was developed using laboratory mice of known chronological age and subsequently applied to wild individuals.

The ecological component includes longitudinal field observations collected during 2023–2024 from two woodland sites near Edinburgh. A factorial field experiment examined the effects of food supplementation and anthelminthic treatment on ecological and ageing-related outcomes.

## Repository structure

```text
├── Data/
│   ├── Processed Data/
│   ├── PanMamManifest_Aligned_ASylvaticus.csv
│   ├── WMClock_3Batches_SampleInfo.csv
│   └── weights_panmamC*.csv
│
├── Scripts/
│   ├── 01_WMClock_Normalization.R
│   ├── 02_WMClock_ClockCreation.R
│   ├── 03_WMClock_EWAS.R
│   └── 04_WMClock_EcologicalAnalyses.R
│
├── Results/
│   ├── 02_WMClock_ClockCreation/
│   ├── 03_WMClock_EWAS/
│   └── 04_WMClock_EcologicalAnalyses/
│
└── README.md
```

### Analysis scripts

| Script | Description |
|---|---|
| `01_WMClock_Normalization.R` | DNA methylation quality control, normalization, batch correction, and annotation |
| `02_WMClock_ClockCreation.R` | Species-specific clock training, validation, pan-mammalian clock comparisons, and tissue analyses |
| `03_WMClock_EWAS.R` | Epigenome-wide association analysis of chronological age and functional enrichment |
| `04_WMClock_EcologicalAnalyses.R` | Ecological correlates of epigenetic age, treatment effects, parasite analyses, and longitudinal epigenetic ageing |

## Reproducing the analyses

Analyses were conducted in R. Scripts are numbered in their intended execution order and use relative file paths.

To reproduce the analyses:

1. Clone or download this repository.
2. Open the associated RStudio project or set the working directory to the repository root.
3. Install the R packages listed at the beginning of each script.
4. Ensure the required input datasets are available in the corresponding `Data/` directories.
5. Run the scripts sequentially.

Scripts automatically create relevant output directories where applicable. Figures and analytical tables are saved within `Results/`, organized by analysis.

Raw methylation IDAT files are required to reproduce the complete normalization workflow. Where these files are unavailable, downstream analyses can be conducted using the processed DNA methylation datasets, provided the required files are accessible.

## Data availability

This repository contains analysis scripts and associated processed data and results, where included.

Availability of raw and processed DNA methylation data is described in the associated manuscript and data repository. Some large input files may be hosted separately from GitHub.

## Citation

If you use these scripts or data, please cite the associated publication:

Wolf et al., 2026 "Epigenetic Clocks Uncover the Natural History of Ageing in Wild Mice" in Molecular Ecology

A permanent archived version of this repository will be available through Zenodo.

## Contact

For questions about the analyses, data, or potential collaborations, please contact the corresponding author through the associated publication.
