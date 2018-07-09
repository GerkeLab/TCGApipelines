match_rules <- tibble::tribble(
   ~liu2018,                              ~cgds,                                  ~merged,
   "bcr_patient_barcode",                 "bcr_patient_barcode",                  "bcr_patient_barcode",
   "type",                                "cancer_type",                          "type",
   "age_at_initial_pathologic_diagnosis", "age",                                  "age_at_initial_pathologic_diagnosis",
   "gender",                              "sex",                                  "gender",
   "race",                                "race",                                 "race",
   "clinical_stage",                      "clin_t_stage",                         "clinical_stage",
   "histological_type",                   "histological_diagnosis",               "histological_type",
   "initial_pathologic_dx_year",          "initial_pathologic_dx_year",           "initial_pathologic_dx_year",
   "birth_days_to",                       "days_to_birth",                        "birth_days_to",
   "vital_status",                        "vital_status",                         "vital_status",
   "tumor_status",                        "tumor_status",                         "tumor_status",
   "last_contact_days_to",                "days_to_last_followup",                "last_contact_days_to",
   "death_days_to",                       "days_to_death",                        "death_days_to",
   "cause_of_death",                      "patient_death_reason",                 "cause_of_death",
   "new_tumor_event_dx_days_to",          "days_to_biochemical_recurrence_first", "new_tumor_event_dx_days_to",
   "treatment_outcome_first_course",      "treatment_outcome_first_course",       "treatment_outcome_first_course",
   "residual_tumor",                      "residual_tumor",                       "residual_tumor",
   "OS",                                  "os_status",                            "OS",
   "OS.time",                             "os_months",                            "OS.time"
)

readr::write_csv(match_rules, here::here("docs/compare-cgds-vs-liu2018/match_rules.csv"))