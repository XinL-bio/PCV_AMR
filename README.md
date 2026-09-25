# Pneumococcal vaccination and antimicrobial resistance

We used ATLAS surveillance data to examine how national introduction of pneumococcal conjugate vaccines (PCVs) relates to antimicrobial non-susceptibility in Streptococcus pneumoniae. The R script produces the study's estimates, tables, and figures.

## Data

The analysis uses the [ATLAS antimicrobial surveillance dataset](https://searchamr.vivli.org/datasetDetails/fromSearch/ce7367e1-c91a-4de1-87e4-6e4d8f834085), available through the [Vivli AMR Register](https://amr.vivli.org/resources/data-request-process-overview/), and a separate country-level PCV introduction-year lookup. Access to ATLAS does not include the PCV lookup. Neither workbook is distributed in this repository; exact results require the same study inputs and country names.

The expected inputs are `atlas_sp.xlsx` (ATLAS isolate records) and `PCV_intro.xlsx`. 

## Run

```bash
Rscript analysis.R
```

## Analysis and outputs

The script estimates changes in non-susceptibility to penicillin, macrolides, and clindamycin, as well as multidrug resistance (MDR), using Sun–Abraham models and a Callaway–Sant’Anna sensitivity analysis. MDR is defined as non-susceptibility to at least three of four antimicrobial classes among isolates tested for at least three classes. Subgroup analyses examine age, specimen type, baseline resistance, country income, and adult PCV programs.

Tables, figures, diagnostic summaries, fitted models, and R session information are written to results/ by default.

## License

The code and documentation are intended for release under the [MIT License](LICENSE). Before publishing, replace the copyright-holder placeholder in `LICENSE` with the appropriate name. The license does not cover ATLAS data, the PCV workbook, or third-party material.
