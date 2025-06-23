---
title: "Mendelian Randomization of CRP and IL-6 in CHD"
author: "Your Name"
date: "2025-06-23"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# 🧬 Tutorial: Investigating CRP and IL-6 in CHD Using Mendelian Randomization

## 📌 Background

Observational studies suggest a causal relationship between elevated CRP levels and coronary heart disease (CHD). However, MR studies using genetic instruments for CRP often fail to support this. IL-6, an upstream regulator of CRP, might be the actual driver. This tutorial demonstrates this with real GWAS data.

## 🧰 1. Setup: Load Packages

```{r packages}
library(TwoSampleMR)
library(ieugwasr)
library(MRInstruments)
library(gwasglue)
library(dplyr)
library(VariantAnnotation)
```

## 📂 2. Accessing GWAS Data

### 🗂 A. Using VCF Files (IEU GWAS local download)

```{r vcf}
vcf_crp <- VariantAnnotation::readVcf("data/ieugwas/CRP_ebi-a-GCST005067_sub.vcf.gz")
vcf_chd <- VariantAnnotation::readVcf("data/ieugwas/CHD_ebi-a-GCST000998_sub.vcf.gz")

crp_exposure <- gwasvcf_to_TwoSampleMR(vcf_crp, "exposure") %>% filter(pval.exposure < 5e-8)
chd_outcome <- gwasvcf_to_TwoSampleMR(vcf_chd, "outcome")

harm_data <- harmonise_data(crp_exposure, chd_outcome)
```

## 🔗 3. Clumping Instruments

```{r clumping}
clump <- ld_clump(
  dplyr::tibble(rsid = harm_data$SNP, pval = harm_data$pval.exposure, id = harm_data$exposure),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  opengwas_jwt = get_opengwas_jwt()
)

harm_clump <- merge(harm_data, clump,
                    by.x = c("SNP", "pval.exposure", "exposure"),
                    by.y = c("rsid", "pval", "id")) %>%
  unique()
```

## 📈 4. MR Analysis: CRP ➡ CHD

```{r mr-analysis}
res <- mr(harm_clump, method_list = c("mr_egger_regression", "mr_ivw"))
print(res)
```

## 🔍 5. Sensitivity Analyses

```{r sensitivity}
mr_heterogeneity(harm_clump)
mr_pleiotropy_test(harm_clump)

# Plots
p1 <- mr_scatter_plot(res, harm_clump)
print(p1)

p2 <- mr_forest_plot(mr_singlesnp(harm_clump))
print(p2)

p3 <- mr_leaveoneout_plot(mr_leaveoneout(harm_clump))
print(p3)

# Report
mr_report(harm_clump, output_path = "figure/")
```

## 🔁 6. Repeat Analysis: IL-6 ➡ CHD

```{r il6-analysis}
exposure <- read_exposure_data(
  filename = "data/gwas_catalogue/IL6-GCST004446.tsv.gz",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
) %>% filter(pval.exposure < 5e-8)
exposure$exposure <- "IL-6"

outcome <- read_outcome_data(
  snps = exposure$SNP,
  filename = "data/gwas_catalogue/CHD_GCST000998.tsv.gz",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)
outcome$outcome <- "CHD"

il6_dat <- harmonise_data(exposure, outcome)

clump_il6 <- ld_clump(
  dplyr::tibble(rsid = il6_dat$SNP, pval = il6_dat$pval.exposure, id = il6_dat$exposure),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  opengwas_jwt = get_opengwas_jwt()
)

harm_clump_il6 <- merge(il6_dat, clump_il6,
                        by.x = c("SNP", "pval.exposure", "exposure"),
                        by.y = c("rsid", "pval", "id")) %>%
  unique()

res_il6 <- mr(harm_clump_il6, method_list = c("mr_egger_regression", "mr_ivw"))
print(res_il6)

mr_heterogeneity(harm_clump_il6)
mr_pleiotropy_test(harm_clump_il6)
mr_report(harm_clump_il6, output_path = "figure/")
```

## 🧠 Summary & Interpretation

| Exposure | Outcome | MR Result | Observational Evidence |
|----------|---------|-----------|--------------------------|
| CRP      | CHD     | No causal effect | Suggests causality |
| IL-6     | CHD     | Causal effect | Consistent with literature |

CRP may not be a causal driver of CHD. IL-6, as an upstream cytokine, likely plays a causal role.

## 📚 References

- Davey Smith & Hemani, 2014. "Mendelian randomization: genetic anchors for causal inference in epidemiological studies." *Hum Mol Genet*.
- Ligthart et al., 2018. *Am J Hum Genet*.
- Georgakis et al., 2020. *Circulation*.
