# Redox / bioenergetics bubbleplot

R script that takes raw biochemical measurements from several species exposed to noise treatments, compares each treatment with its controls, and summarises all comparisons in bubbleplots (assay × species, split by time point and treatment).

## Usage

```bash
Rscript redox_bubbleplot.R path/to/input.xlsx path/to/output_folder
```

When you run it interactively (for example in RStudio), the script asks for the input file and the output folder. Missing R packages (`tidyverse`, `readxl`, `emmeans`, `broom`, `scales`) are installed on the first run.

Edit the **Settings** block at the top of the script to change treatment codes, species and assay labels, time points, figure size, or to pool life stages (see below).

## Analysis

1. **Technical replicates.** All rows that share `sample_id`, `species`, `treatment`, `time_h`, `assay` and `batch` are averaged, giving one value per biological sample. Statistics are always computed on biological samples, never on technical replicates.
2. **P/O ratio.** OxPhos efficiency is computed for each biological sample as ATP synthesis / OCR. 
3. **Normalisation.** Each sample is divided by the mean of the controls of the same species, assay, time point and batch. If a batch has no controls, the mean of all controls of that species, assay and time point is used instead. The result is the fold change (FC), then log2-transformed (log2FC).
4. **Statistics.** For each species × assay × time point, a one-way linear model (`log2FC ~ treatment`) is fitted, followed by treatment-vs-control contrasts (*emmeans*, `trt.vs.ctrl`). These are t-tests on the residual variance pooled across all groups, unadjusted at this stage.
5. **Multiple testing.** p-values are adjusted with the Benjamini–Hochberg procedure (FDR) across all tests of a figure set (see *Outputs*). Because the correction is done separately for each set, the same comparison can have slightly different adjusted p-values in the two sets.
6. **Effect size.** Cohen's d and Hedges' g are computed on log2FC for each treatment vs its controls: difference between means divided by the pooled SD of the two groups. Hedges' g applies the small-sample correction J = 1 − 3 / (4·df − 1), with df = n₁ + n₂ − 2.

### Life stages

Rows whose `notes` contain "larva" (case-insensitive) are analysed as a **separate unit** from the adults/juveniles of the same species. They have their own controls, tests and row in the figure, even when the `species` code is the same. Setting `POOL_LIFE_STAGES <- TRUE` analyses all stages of a species together; this is intended as a proof-of-concept comparison only. Normalisation is still batch-matched, so each sample is compared with controls of its own batch.

## Outputs

Written to the chosen output folder:

| File | Content |
|---|---|
| `01_PO_calculated.csv` | P/O ratio for each biological sample |
| `02_data_normalized_log2FC.csv` | Biological-sample values with control reference, FC and log2FC |
| `03_stats_FDR_effectsizes_<set>.csv` | One row per comparison: estimate (difference in log2FC), SE, df, t, p, FDR-adjusted p, Cohen's d, Hedges' g, n |
| `figures/bubbleplot_log2FC_padj_<set>.pdf/.png` | Colour = mean log2FC; size = FDR-adjusted p-value class |
| `figures/bubbleplot_hedgesG_padj_cap_<set>.pdf/.png` | Colour = Hedges' g, scale capped at the 90th percentile of \|g\| |
| `figures/bubbleplot_hedgesG_padj_nocap_<set>.pdf/.png` | Colour = Hedges' g, scale set by the largest \|g\| |

`<set>` is `noLarvae` (adults/juveniles only) and, if the input contains larvae, `withLarvae` (all units). With `POOL_LIFE_STAGES <- TRUE` there is one set, `pooledStages`.

In every figure, bubbles with FDR-adjusted p > 0.05 have a grey outline, and the others a black outline.

## Input format

One row per measurement (long format), as `.xlsx` (first sheet) or `.csv`. `input_template.csv` contains the header row to fill in.

| Column | Description | Example |
|---|---|---|
| `sample_id` | Biological sample (animal). Rows with the same `sample_id` in the same species, treatment, time, assay and batch are technical replicates | `CL_13` |
| `species` | Species code; must be listed in `SPECIES_LABELS` | `Clepadiformis` |
| `treatment` | `Control` or a treatment code listed in `TREATMENT_RENAME` | `MSFD-S` |
| `time_h` | Time point in hours; must be one of `TIME_POINTS_H` (default 1, 24) | `24` |
| `assay` | Assay name; must be one of the names in `ASSAY_LABELS` except `P/O`. P/O needs both `ATP synthesis` and `OCR` | `CAT` |
| `value` | Measured value (numeric, `.` as decimal separator) | `229.14` |
| `unit` | Unit of `value`, for documentation only | `mU/mg` |
| `replicate_id` | Technical replicate label, for documentation only | `R1` |
| `batch` | Experimental batch. Each batch should include its own controls | `CLAV_1` |
| `notes` | Free text; "Larva" marks larval samples | `Whole adult` |

Before computing anything, the script checks that all columns are present, that `time_h` and `treatment` contain only allowed values, and that every species is listed in the settings. If a check fails, it stops and lists the offending rows.

Each species × assay × time point needs controls and at least two biological samples per group for the tests and effect sizes.

## Notes

- Hedges' g is the difference divided by the variability among animals. With few animals per group (n = 3) and very homogeneous samples, g can become very large even when the change is modest.
- PDF figures are written with cairo when it is available, with Quartz on macOS when cairo is not, and with the base PDF device otherwise.
