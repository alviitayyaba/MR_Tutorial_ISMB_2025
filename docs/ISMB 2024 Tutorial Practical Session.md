# ISMB 2024 Tutorial Practical Session: Mendelian Randomization (MR)

---

## 📅 Agenda (Hands-On Tutorial: 10:15–13:00)

### 🔧 10:15–10:45 — Setup and Getting Started

* **Goals**:

  * Introduce the agenda and expected outcomes.
  * Ensure all participants are set up with R, RStudio, and required packages.
* **Activities**:

  * Guide installation of R packages: `TwoSampleMR`, `ieugwasr`, `tidyverse`, `MRInstruments`.
  * Explore public GWAS data sources:

    * [GWAS Catalog](https://www.ebi.ac.uk/gwas/)
    * [IEU Open GWAS Project](https://gwas.mrcieu.ac.uk/)
  * I will demonstrate access and retrieval of GWAS summary statistics from two sources.
    * GWAS catalogue 
    * ieugwasr databse


### 🧹 11:00-11:15 — Data Harmonization and Instrument Selection

* **Goals**:

  * Prepare exposure and outcome datasets for MR.
  * Understand the importance of LD clumping and harmonization.
* **Activities**:

  * Perform LD clumping to select independent variants.
  * Harmonize allele coding between datasets.
  * Subset GWAS by:

    * Biological significance (e.g. selecting SNPs near IL6R)
    * Statistical significance (p-value threshold)
* **Challenges Addressed**:

  * LD, Strand mismatches, palindromic SNPs, allele frequency conflicts.
* **Outcome**:

  * Clean, harmonized datasets ready for MR analysis.

### 🧪 11:15–12:30 — Mendelian Randomization Analyses

* **Goals**:

  * Conduct MR and apply robust methods to address violations of assumptions.
* **Case Study 1**: CRP ➝ Coronary Heart Disease

  * Purpose: Show how MR reveals a *non-causal* relationship.
  * Demonstrates how observational correlations can be due to confounding.
* **Case Study 2**: IL6 ➝ Coronary Heart Disease

  * Purpose: Show how MR identifies a *causal* effect using biologically informed instruments.
  * Demonstrates the value of using cis-variants and known biology.
* **Activities**:
Taking both of these case studies. I will demonstrate the application of MR different methods within TwoSampleMR pacakge. 
  * Primary MR (Inverse Variance Weighted method)
  * Robust MR methods:

    * MR-Egger
    * Weighted Median
    * Mode-based estimate (optional)
  * Visualization and interpretation of findings:

    * Forest plots
    * Scatter plots
    * Funnel plots
  * Assess heterogeneity and pleiotropy:

    * Cochran's Q, Egger intercept, MR-PRESSO (if time permits)

### 🔍 12:30–13:00 — Detecting and Handling Violations

* **Goals**:

  * Identify invalid instruments and their effects on results.
  * Demonstrate iterative re-analysis after correction.
* **Activities**:

  * Identify outlier SNPs and perform leave-one-out analysis.
  * Re-run MR after removing invalid IVs.
  * Discuss trade-offs: sensitivity vs. specificity.
* **Challenges Addressed**:

  * Pleiotropy, weak instrument bias, invalid IVs.
* **Outcome**:

  * Learn practical strategies to deal with common pitfalls in MR studies.

---

## 🧠 Learning Objectives Summary

By the end of the session, participants will be able to:

* Access and prepare GWAS summary statistics for MR.
* Perform basic and advanced MR analyses in R.
* Detect and address issues like pleiotropy, heterogeneity, and invalid IVs.
* Interpret MR results critically and in context.

---


## Suggested Structure for Theory Component

Based on this practical outline, I recommend the theory sessions focus on building a strong conceptual foundation and awareness of challenges in MR design. The ideal breakdown:

###  First 30 minutes:

* Introduction to Mendelian Randomization:

  * Core idea: Using genetic variants as instruments - analogy to a RCT. 
  * Classic examples where MR contradicted observational findings (e.g. CRP and CHD).
* Assumptions of MR:
  * Relevance, independence, exclusion restriction.
* Overview of GWAS data and one-sample vs. two-sample MR.

###  Next 45 minutes:

* Common pitfalls and challenges:

  * Confounding due to pleiotropy.
  * Horizontal vs. vertical pleiotropy.
  * Weak instruments.
* Methodological responses:

  * MR-Egger regression: intuition and mathematical formulation.
  * Heterogeneity testing (e.g., Cochran’s Q).
  





