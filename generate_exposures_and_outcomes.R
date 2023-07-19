pacman::p_load(tidyverse,vroom,data.table,TwoSampleMR)

il6_path_exposures <- vroom("https://raw.githubusercontent.com/gushamilton/il6-sepsis/main/data/harmonised_data_final.tsv")%>%
  mutate(beta.exposure = if_else(SNP == "rs3093077", -beta.exposure, beta.exposure)) %>%
  mutate(beta.exposure = if_else(SNP == "rs1800947", -beta.exposure, beta.exposure)) %>%
  filter(exposure == "cisCRP" | exposure == "cisIL6R") %>%
  select(SNP, contains("exposure")) %>%
  mutate(beta.exposure = -beta.exposure) %>%
  distinct()

IL6R_start <- 154377819 - 3e5
IL6R_end <- 154441926 +3e5
IL6R_end
tnf_start <- 6437923-3e5 
tnf_end <- 6451280 + 3e5
il1ra_start <- 113875548
il1ra_end <- 113891591
crp  <- data.table::fread("~/data/CRP/gwas/GCST90029070_buildGRCh37.tsv.gz") 

tnf_exposures <- crp %>%
  filter(chromosome == 12 & (base_pair_location > tnf_start) & (base_pair_location < tnf_end)) %>%
  select(SNP = variant_id, effect_allele.exposure = effect_allele, other_allele.exposure = other_allele, beta.exposure = beta, se.exposure = standard_error, pval.exposure = p_value) %>%
  as_tibble() %>%
  mutate(eaf.exposure = NA) %>%
  mutate(exposure ="TNF", id.exposure = "TNF") 

tnf_exposures_clumped <- tnf_exposures%>%
  mutate(pval.exposure = as.numeric(pval.exposure)) %>%
   filter(pval.exposure < 5e-8)  %>%
  clump_data(clump_r2 = 0.5)


il1ra <- crp %>%
  filter(chromosome == 2 & (base_pair_location > il1ra_start) & (base_pair_location < il1ra_end)) %>%
  select(SNP = variant_id, effect_allele.exposure = effect_allele, other_allele.exposure = other_allele, beta.exposure = beta, se.exposure = standard_error, pval.exposure = p_value) %>%
  as_tibble() %>%
  mutate(eaf.exposure = NA) %>%
  mutate(exposure ="IL1RA", id.exposure = "Il1RA") 

il1ra_clumped <- il1ra %>%
  mutate(pval.exposure = as.numeric(pval.exposure)) %>%
  filter(pval.exposure < 5e-8)  %>%
  clump_data(clump_r2 = 0.2)

all_exposures <- bind_rows(tnf_exposures_clumped, il6_path_exposures, il1ra_clumped)
all_exposures %>%
  count(exposure)

ao <-available_outcomes()

lipid_concs <- ao %>%
  filter(str_detect(id, "met-d")) %>%
  filter(str_detect(trait, "Conce|Apo")) %>%
  pull(id)
outcomes_met <- extract_outcome_data(all_exposures$SNP, lipid_concs)

# GLCC


# # use tabix
glcc <- vroom("glcc/glcc_output", col_names = F)
glcc_outcomes <- glcc %>%
  transmute(SNP = X1,
            chrom.outcome = X2,
            pos.outcome = X3,
            other_allele.oucome = X4,
            effect_allele.outcome = X5,
            N = X6,
            eaf.outcome = X8,
            beta.outcome = X9,
            se.outcome = X10,
            pval.outcome = X12,
            outcome = X15
  ) %>%
  mutate(outcome = paste(str_remove(outcome, "_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results"), "GLCC")) %>%
  mutate(id.outcome = outcome)

outcomes <- bind_rows(glcc_outcomes, outcomes_met)
            
dat <- harmonise_data(all_exposures, outcomes, action = 1)
write_tsv(dat, "harmonised_exposure_outcome.dat")
