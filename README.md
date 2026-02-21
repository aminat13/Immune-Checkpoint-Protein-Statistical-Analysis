# Immune Checkpoint Protein Statistical Analysis

This repository contains coursework-based statistical analysis of immune checkpoint protein expression data in ovarian cancer patients treated with Rucaparib or Niraparib.

The analyses were conducted as part of a biostatistics module and were structured around specific statistical questions. The assignments focused on parametric testing, non-parametric testing, and effect size/power analysis.

The code reflects the original analytical approach used to answer the assignment briefs.

---

# Dataset Overview

The dataset contains expression levels of four immune checkpoint proteins:

- PD-1  
- PDL-1  
- CTLA-4  
- SIRP  

Measured across treatment arms (Rucaparib and Niraparib), across timepoints (pre-, during-, post-treatment), and alongside clinical variables including ER status, treatment response, and BMI.

---

# Assignment 1 – Parametric Hypothesis Testing

## Questions Addressed

1. Is there a statistically significant difference in the expression levels of the four immune checkpoint proteins between treatment arms before treatment (analyze each marker separately)?

2. Are oestrogen receptor (ER) positive samples more likely to be non-responders?

3. Is there a statistically significant difference between the expression levels of the four immune checkpoint proteins within the Rucaparib treatment arm (pre-treatment), and where does that difference lie?

## Statistical Methods Used

- Shapiro-Wilk test (normality assessment)
- Levene’s test (homogeneity of variance)
- Independent Student’s t-tests
- Wilcoxon rank-sum test (where assumptions were not met)
- Fisher’s Exact test
- Kruskal-Wallis test
- Post-hoc analysis with Benjamini-Hochberg correction

Script: `scripts/parametric_hypothesis_testing.R`

---

# Assignment 2 – Non-Parametric Testing

## Questions Addressed

1. Is there a statistically significant difference in protein expression levels pre-treatment and post-treatment (analyzing each marker separately)?

2. Is there a significant increase in the number of underweight patients post-treatment (BMI < 20)?

3. Is there a statistically significant difference in protein expression across the course of the trial (pre-, during-, and post-treatment), and where does that difference lie?

## Statistical Methods Used

- Wilcoxon signed-rank test
- Friedman test for repeated measures
- Post-hoc pairwise comparisons
- McNemar’s test for paired categorical data

Script: `scripts/non_parametric_testing.R`

---

# Assignment 3 – Effect Size and Power Analysis

## Questions Addressed

1. For markers that were not statistically significant (pre-treatment), how many samples would be required to achieve 80% power to detect a statistically significant difference?

2. Are the immune checkpoint markers affected differently by anti-PARP treatments (post-treatment cohort)? Where differences were not significant, how many samples would be required to detect an effect with 80% power?

3. Is there a significant increase in the number of underweight patients post-treatment? If not, how many samples would be required to achieve 80% power?

4. Is obesity prevalence in the ovarian cancer cohort significantly different from the Irish adult population (23%)? What sample size would be required to detect this difference with 80% power?

## Statistical Methods Used

- Cohen’s d
- Cohen’s f
- Cohen’s h
- Cohen’s w
- Power calculations using the `pwr` package
- ANOVA and Tukey post-hoc testing (where applicable)

Script: `scripts/effect_size_and_power_analysis.R`

