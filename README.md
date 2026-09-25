# Pneumococcal vaccination and antimicrobial resistance

This repository contains the R analysis for a study of national pneumococcal conjugate vaccine (PCV) introduction and antimicrobial non-susceptibility in *Streptococcus pneumoniae*. It generates the study's estimates, tables, and figures from isolate-level surveillance data and country-specific PCV introduction years.

## Data

The analysis uses the [ATLAS antimicrobial surveillance dataset](https://searchamr.vivli.org/datasetDetails/fromSearch/ce7367e1-c91a-4de1-87e4-6e4d8f834085), available through the [Vivli AMR Register](https://amr.vivli.org/resources/data-request-process-overview/), and a separate country-level PCV introduction-year lookup. Access to ATLAS does not include the PCV lookup. Neither workbook is distributed in this repository; exact results require the same study inputs and country names.

The expected inputs are `atlas_sp.xlsx` (ATLAS isolate records) and `PCV_intro.xlsx` (columns `Country` and `Year_all`). The script checks the necessary ATLAS columns and the study cohort before analysis.

## Run

Use R with `readxl`, `dplyr`, `ggplot2`, `fixest`, `did`, `tibble`, and `patchwork`. The analysis was checked with R 4.4.2, `fixest` 0.14.2, and `did` 2.5.1.

With both workbooks in the working directory, run:

```bash
Rscript analysis.R
```

To keep the workbooks elsewhere, set `ATLAS_SP_XLSX` and `PCV_INTRO_XLSX` to their full paths before running. Set `PCV_OUTPUT_DIR` to change the default `results/` output directory.

## Analysis and outputs

The script estimates changes in penicillin, macrolide, clindamycin, and multidrug non-susceptibility using Sun–Abraham models, with a Callaway–Sant’Anna sensitivity analysis. Multidrug resistance (MDR) is defined as non-susceptibility to at least three of four antimicrobial classes among isolates tested for at least three classes. Subgroup analyses examine age, specimen type, baseline resistance, country income, and adult PCV programs.

Tables, figures, diagnostic summaries, fitted models, and R session information are written to `results/` by default. Figure 3 reports effects as proportions (for example, −0.5 equals −50 percentage points). Do not publish derived country-level outputs without checking the applicable data-use terms.

## License

The code and documentation are intended for release under the [MIT License](LICENSE). Before publishing, replace the copyright-holder placeholder in `LICENSE` with the appropriate name. The license does not cover ATLAS data, the PCV workbook, or third-party material.
