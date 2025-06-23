library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(VariantAnnotation)
devtools::install_github("mrcieu/gwasglue")
library(gwasglue)

### read data with TwosampleMR pacakge from ieugwasr databses reading from local download
vcf_crp <- VariantAnnotation::readVcf("data/ieugwas/CRP_ebi-a-GCST005067_sub.vcf.gz")
crp_exposure <- gwasvcf_to_TwoSampleMR(vcf_crp, "exposure")
crp_exposure <- crp_exposure %>% filter(pval.exposure < 5*10^-8)
vcf_chd <- VariantAnnotation::readVcf("data/ieugwas/CHD_ebi-a-GCST000998_sub.vcf.gz")
chd_outcome <- gwasvcf_to_TwoSampleMR(vcf_chd, "outcome")
harm_data <- harmonise_data(exposure_dat = crp_exposure, outcome_dat = chd_outcome)

#clump <- ld_clump(
#  dplyr::tibble(rsid=harm_data$SNP, pval=harm_data$pval.exposure, id=harm_data$exposure),
#  pop = "EUR",
#  plink_bin = genetics.binaRies::get_plink_binary(),
#  bfile = "/wins/donertas/personal/Tayyaba/1kg.v3/EUR"
#)


### read data with MRInstruments pacakge from Mr base databses using API
library(MRInstruments)
vignette("MRBase")
data(gwas_catalog)
exposure_gwas <- subset(gwas_catalog, grepl("GCST007615", STUDY.ACCESSION) &
                          Phenotype_simple == "C-reactive protein levels")
exposure_gwas<-exposure_gwas[exposure_gwas$pval<5*10^-8,]
exposure_data<-format_data(exposure_gwas)
head(exposure_data)

ao <- available_outcomes()
head(ao)
##according to new API for Open-GWAS this rquires limit on account usage and usually does not work with trail accounts after the limit is spent
chd_outcome <- extract_outcome_data(
  snps = exposure_data$SNP,
  outcomes = 'ieu-a-7'
)
### read data with TwosampleMr from gwas catalogue databses from local download 
exposure <- read_exposure_data(
    filename = "data/gwas_catalogue/CRP-GCST005067_sub.tsv.gz",
    sep = "\t",
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )

exposure$exposure <- "CRP"
exposure <- exposure %>%
  filter(pval.exposure < 5*10^-8)
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
dat <- harmonise_data(
  exposure_dat =  exposure,
  outcome_dat = outcome
)
clump <- ld_clump(
  dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$exposure),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  opengwas_jwt = get_opengwas_jwt(),
  bfile = NULL,
  plink_bin = NULL
)
harm_clump <- merge(dat, clump, by.x = c("SNP", "pval.exposure", "exposure"), by.y = c("rsid", "pval", "id")) %>% unique()
res <- mr(harm_clump, method_list = c("mr_egger_regression", "mr_ivw"))
res
mr_heterogeneity(harm_clump)
mr_pleiotropy_test(harm_clump)
p1 <- mr_scatter_plot(res, harm_clump)
p1
res_single <- mr_singlesnp(harm_clump)
p2 <- mr_forest_plot(res_single)
p2
res_loo <- mr_leaveoneout(harm_clump)
p3 <- mr_leaveoneout_plot(res_loo)
p3
mr_report(harm_clump, output_path = "figure/")

outcome <- read_outcome_data(
  snps = exposure$SNP,
  filename = "data/gwas_catalogue/IL6-GCST004446.tsv.gz",
  sep = "\t",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)
outcome$outcome <- "IL-6"
dat <- harmonise_data(
  exposure_dat =  exposure,
  outcome_dat = outcome
)
clump <- ld_clump(
  dplyr::tibble(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$exposure),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  opengwas_jwt = get_opengwas_jwt(),
  bfile = NULL,
  plink_bin = NULL
)
harm_clump <- merge(dat, clump, by.x = c("SNP", "pval.exposure", "exposure"), by.y = c("rsid", "pval", "id")) %>% unique()
res <- mr(harm_clump, method_list = c("mr_egger_regression", "mr_ivw"))
res
mr_heterogeneity(harm_clump)
mr_pleiotropy_test(harm_clump)
p1 <- mr_scatter_plot(res, harm_clump)
p1
res_single <- mr_singlesnp(harm_clump)
p2 <- mr_forest_plot(res_single)
p2
res_loo <- mr_leaveoneout(harm_clump)
p3 <- mr_leaveoneout_plot(res_loo)
p3
mr_report(harm_clump, output_path = "figure/")



exposure_gwas <- subset(gwas_catalog, grepl("GCST004446", STUDY.ACCESSION) &
                          Phenotype == "Interleukin-6 levels")
exposure_gwas <- exposure_gwas[exposure_gwas$pval < 5*10^-8, ]
