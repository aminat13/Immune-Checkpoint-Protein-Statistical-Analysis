# Immune Checkpoint Protein Statistical Analysis in R

This project performs a comprehensive statistical analysis of PD-1, PDL-1, CTLA-4, and SIRP
expression in ovarian cancer patients treated with Rucaparib or Niraparib.

## Methods Used
- Exploratory data analysis (histograms, boxplots)
- Assumption testing (Shapiro-Wilk, Levene’s test)
- Parametric & non-parametric testing
- Multiple comparison correction (Benjamini-Hochberg)
- Longitudinal analysis (Friedman test)
- Categorical analysis (Fisher’s exact, McNemar’s test)
- Effect size estimation (Cohen’s d, f, h, w)
- Power analysis using `pwr`

## Key Findings
- PD-1 expression differs between treatment arms
- Significant longitudinal changes observed across timepoints
- ER status associated with treatment response
- Power calculations indicate oversampling in some comparisons
