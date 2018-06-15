# ---- Library ----
library(tidyverse)
if (!requireNamespace('pacman', quietly = TRUE)) install.packages("pacman")
pacman::p_load("cgdsr")

# ---- Setup ----
short_name <- "tcga-cgds_prad"
cgds <- CGDS("http://www.cbioportal.org/public-portal/")

# ---- Download Data ----
prad_tcga <- as_tibble(getCancerStudies(cgds)) %>% 
   filter(str_detect(name, "Prostate")) %>% 
   filter(str_detect(name, "TCGA, Provisional"))

prad_tcga_caselist <- getCaseLists(cgds, prad_tcga$cancer_study_id) %>% 
   as_tibble()

prad_tcga_clinicaldata <- prad_tcga_caselist %>% 
   filter(str_detect(case_list_id, "_all$")) %>% 
   pull(case_list_id) %>% 
   .[1] %>% 
   getClinicalData(cgds, .) %>% 
   tibble::rownames_to_column("id") %>% 
   as_tibble() %>% 
   janitor::clean_names()

saveRDS(prad_tcga_clinicaldata, here::here("data", paste0(short_name, ".rds")))
readr::write_csv(prad_tcga_clinicaldata, here::here("data", paste0(short_name, ".csv")))
