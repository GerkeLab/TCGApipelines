
## Clinical data exploration example: prostate

2018-06-14

-----

Data used here was downloaded from TCGA’s pan-cancer clinical data
resource from Supplementary Table 1 of [this
paper](https://www.cell.com/cell/fulltext/S0092-8674\(18\)30229-0).

``` r
library(here)
library(tidyverse)

load(here("data", "tcga-liu2018_clinData.RData"))

# filter to PRAD
dat <- dat %>% filter(type=="PRAD")
glimpse(dat)
```

    ## Observations: 500
    ## Variables: 34
    ## $ X__1                                <chr> "7867", "7868", "7869", "7...
    ## $ bcr_patient_barcode                 <chr> "TCGA-2A-A8VL", "TCGA-2A-A...
    ## $ type                                <chr> "PRAD", "PRAD", "PRAD", "P...
    ## $ age_at_initial_pathologic_diagnosis <dbl> 51, 57, 47, 52, 70, 54, 69...
    ## $ gender                              <chr> "MALE", "MALE", "MALE", "M...
    ## $ race                                <chr> "[Not Available]", "[Not A...
    ## $ ajcc_pathologic_tumor_stage         <chr> "[Not Applicable]", "[Not ...
    ## $ clinical_stage                      <chr> "[Not Applicable]", "[Not ...
    ## $ histological_type                   <chr> "Prostate Adenocarcinoma A...
    ## $ histological_grade                  <chr> "[Not Available]", "[Not A...
    ## $ initial_pathologic_dx_year          <dbl> 2010, 2010, 2011, 2010, 20...
    ## $ menopause_status                    <chr> "[Not Available]", "[Not A...
    ## $ birth_days_to                       <dbl> -18658, -20958, -17365, -1...
    ## $ vital_status                        <chr> "Alive", "Alive", "Alive",...
    ## $ tumor_status                        <chr> "TUMOR FREE", "TUMOR FREE"...
    ## $ last_contact_days_to                <dbl> 621, 1701, 1373, 671, 1378...
    ## $ death_days_to                       <dbl> NA, NA, NA, NA, NA, NA, NA...
    ## $ cause_of_death                      <chr> "[Not Available]", "[Not A...
    ## $ new_tumor_event_type                <chr> NA, NA, NA, NA, NA, NA, "B...
    ## $ new_tumor_event_site                <chr> NA, NA, NA, NA, NA, NA, NA...
    ## $ new_tumor_event_site_other          <chr> NA, NA, NA, NA, NA, NA, NA...
    ## $ new_tumor_event_dx_days_to          <dbl> NA, NA, NA, NA, NA, NA, 19...
    ## $ treatment_outcome_first_course      <chr> "[Not Available]", "[Not A...
    ## $ margin_status                       <chr> NA, NA, NA, NA, NA, NA, NA...
    ## $ residual_tumor                      <chr> NA, NA, NA, NA, NA, NA, NA...
    ## $ OS                                  <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    ## $ OS.time                             <dbl> 621, 1701, 1373, 671, 1378...
    ## $ DSS                                 <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0,...
    ## $ DSS.time                            <dbl> 621, 1701, 1373, 671, 1378...
    ## $ DFI                                 <dbl> NA, NA, NA, NA, NA, NA, NA...
    ## $ DFI.time                            <dbl> NA, NA, NA, NA, NA, NA, NA...
    ## $ PFI                                 <dbl> 0, 0, 0, 0, 0, 0, 1, 0, 0,...
    ## $ PFI.time                            <dbl> 621, 1701, 1373, 671, 1378...
    ## $ Redaction                           <chr> NA, NA, NA, NA, NA, NA, NA...

Try to figure out how time to event variables are created.

``` r
# OS.time and DSS.time are the same thing
all.equal(dat$OS.time, dat$DSS.time)
```

    ## [1] TRUE

``` r
# new_tumor_event_type seems to denote what sort of event was recorded
table(dat$new_tumor_event_type, useNA="ifany")
```

    ## 
    ## Biochemical evidence of disease              Distant Metastasis 
    ##                              70                               6 
    ##         Locoregional Recurrence               New Primary Tumor 
    ##                               8                               7 
    ##                            <NA> 
    ##                             409

``` r
# 5 patients with unequal DFI and PFI times when DFI is not mising, all have new primaries
dat %>% 
   filter(DFI.time != PFI.time) %>% 
   select(bcr_patient_barcode, new_tumor_event_type, DFI.time, PFI.time) %>% 
   knitr::kable() 
```

| bcr\_patient\_barcode | new\_tumor\_event\_type | DFI.time | PFI.time |
| :-------------------- | :---------------------- | -------: | -------: |
| TCGA-CH-5791          | New Primary Tumor       |     1004 |      396 |
| TCGA-HC-7080          | New Primary Tumor       |     1106 |      894 |
| TCGA-J4-A83N          | New Primary Tumor       |      992 |      433 |
| TCGA-VP-A87B          | New Primary Tumor       |     2722 |     2473 |
| TCGA-XK-AAIW          | New Primary Tumor       |     1218 |      427 |

``` r
# what about those with BCR?
dat %>% 
   filter(new_tumor_event_type=="Biochemical evidence of disease") %>% 
   select(bcr_patient_barcode, new_tumor_event_type, DFI.time, PFI.time) %>% 
   knitr::kable() 
```

| bcr\_patient\_barcode | new\_tumor\_event\_type         | DFI.time | PFI.time |
| :-------------------- | :------------------------------ | -------: | -------: |
| TCGA-2A-A8W3          | Biochemical evidence of disease |       NA |      198 |
| TCGA-EJ-8472          | Biochemical evidence of disease |       NA |      196 |
| TCGA-EJ-A46F          | Biochemical evidence of disease |       NA |      215 |
| TCGA-EJ-A65F          | Biochemical evidence of disease |       NA |       75 |
| TCGA-EJ-A7NN          | Biochemical evidence of disease |      197 |      197 |
| TCGA-EJ-A8FP          | Biochemical evidence of disease |      117 |      117 |
| TCGA-EJ-A8FS          | Biochemical evidence of disease |      216 |      216 |
| TCGA-G9-6332          | Biochemical evidence of disease |     1180 |     1180 |
| TCGA-G9-6339          | Biochemical evidence of disease |     1634 |     1634 |
| TCGA-G9-6498          | Biochemical evidence of disease |     1342 |     1342 |
| TCGA-HC-7079          | Biochemical evidence of disease |       NA |      380 |
| TCGA-HC-7213          | Biochemical evidence of disease |       NA |      170 |
| TCGA-HC-7232          | Biochemical evidence of disease |      766 |      766 |
| TCGA-HC-7738          | Biochemical evidence of disease |       NA |      420 |
| TCGA-HC-A9TE          | Biochemical evidence of disease |       NA |      216 |
| TCGA-HC-A9TH          | Biochemical evidence of disease |       NA |      351 |
| TCGA-HI-7168          | Biochemical evidence of disease |       NA |     2505 |
| TCGA-J4-A67N          | Biochemical evidence of disease |       NA |      442 |
| TCGA-J4-A67S          | Biochemical evidence of disease |      708 |      708 |
| TCGA-J4-A6G3          | Biochemical evidence of disease |       NA |      618 |
| TCGA-J4-A83M          | Biochemical evidence of disease |      521 |      521 |
| TCGA-J4-AATZ          | Biochemical evidence of disease |       79 |       79 |
| TCGA-J9-A52B          | Biochemical evidence of disease |       57 |       57 |
| TCGA-J9-A8CL          | Biochemical evidence of disease |       NA |      132 |
| TCGA-J9-A8CM          | Biochemical evidence of disease |       NA |      344 |
| TCGA-KC-A4BL          | Biochemical evidence of disease |       NA |      193 |
| TCGA-KC-A4BR          | Biochemical evidence of disease |     1022 |     1022 |
| TCGA-KC-A4BV          | Biochemical evidence of disease |     1328 |     1328 |
| TCGA-KK-A6E0          | Biochemical evidence of disease |      941 |      941 |
| TCGA-KK-A6E7          | Biochemical evidence of disease |       NA |      925 |
| TCGA-KK-A7AQ          | Biochemical evidence of disease |     1217 |     1217 |
| TCGA-KK-A7AU          | Biochemical evidence of disease |       NA |      207 |
| TCGA-KK-A7AY          | Biochemical evidence of disease |     1124 |     1124 |
| TCGA-KK-A7B0          | Biochemical evidence of disease |       NA |      606 |
| TCGA-KK-A7B2          | Biochemical evidence of disease |       NA |      692 |
| TCGA-KK-A7B3          | Biochemical evidence of disease |      294 |      294 |
| TCGA-KK-A7B4          | Biochemical evidence of disease |       NA |      637 |
| TCGA-KK-A8I7          | Biochemical evidence of disease |       NA |     1088 |
| TCGA-KK-A8I9          | Biochemical evidence of disease |      940 |      940 |
| TCGA-KK-A8IC          | Biochemical evidence of disease |       NA |     1060 |
| TCGA-KK-A8IF          | Biochemical evidence of disease |       NA |      648 |
| TCGA-KK-A8II          | Biochemical evidence of disease |       NA |      626 |
| TCGA-M7-A722          | Biochemical evidence of disease |       NA |      559 |
| TCGA-V1-A8MM          | Biochemical evidence of disease |      990 |      990 |
| TCGA-V1-A9O5          | Biochemical evidence of disease |      124 |      124 |
| TCGA-V1-A9O7          | Biochemical evidence of disease |       NA |      819 |
| TCGA-V1-A9OL          | Biochemical evidence of disease |      105 |      105 |
| TCGA-V1-A9OT          | Biochemical evidence of disease |       NA |      292 |
| TCGA-VN-A88R          | Biochemical evidence of disease |       NA |      512 |
| TCGA-VP-A878          | Biochemical evidence of disease |       NA |       98 |
| TCGA-VP-A87D          | Biochemical evidence of disease |     1194 |     1194 |
| TCGA-VP-A87K          | Biochemical evidence of disease |      533 |      533 |
| TCGA-XK-AAJR          | Biochemical evidence of disease |      131 |      131 |
| TCGA-YL-A8HK          | Biochemical evidence of disease |     1376 |     1376 |
| TCGA-YL-A8HM          | Biochemical evidence of disease |     1473 |     1473 |
| TCGA-YL-A8HO          | Biochemical evidence of disease |       NA |     1068 |
| TCGA-YL-A8S8          | Biochemical evidence of disease |       NA |      679 |
| TCGA-YL-A8SB          | Biochemical evidence of disease |       NA |     1384 |
| TCGA-YL-A8SC          | Biochemical evidence of disease |       NA |      152 |
| TCGA-YL-A8SI          | Biochemical evidence of disease |     1423 |     1423 |
| TCGA-YL-A8SJ          | Biochemical evidence of disease |      752 |      752 |
| TCGA-YL-A8SP          | Biochemical evidence of disease |       NA |     2036 |
| TCGA-YL-A8SQ          | Biochemical evidence of disease |       NA |      329 |
| TCGA-YL-A9WK          | Biochemical evidence of disease |       NA |     1009 |
| TCGA-YL-A9WL          | Biochemical evidence of disease |      740 |      740 |
| TCGA-YL-A9WX          | Biochemical evidence of disease |       NA |     1506 |
| TCGA-YL-A9WY          | Biochemical evidence of disease |       NA |      765 |
| TCGA-ZG-A9L2          | Biochemical evidence of disease |       NA |      180 |
| TCGA-ZG-A9L6          | Biochemical evidence of disease |       NA |      664 |
| TCGA-ZG-A9L9          | Biochemical evidence of disease |       NA |       51 |

``` r
# look at time to BCR in days and then in months
BCRdat <- dat %>% filter(new_tumor_event_type=="Biochemical evidence of disease")
summary(BCRdat$PFI.time)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    51.0   209.0   622.0   679.6  1018.8  2505.0

``` r
summary(BCRdat$PFI.time/30.5)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   1.672   6.852  20.393  22.282  33.402  82.131
