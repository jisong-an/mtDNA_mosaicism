################################
#
#   Shared variant classification
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

library(tidyverse)


## 1. make dataframe of shared variants ---------------------------------

df_shared_var <- df %>% 
  filter(type=="clonal_normal") %>%
  select(project,patient,sample,VAR) %>% 
  unique() %>% 
  filter(!is.na(VAR)) %>% 
  mutate(total_count=nrow(sample_info_final_v3 %>% filter(type=="clonal_normal"))) %>% 
  add_count(VAR, name="VAR_total_count")  %>% 
  merge(sample_info_final_v3 %>% filter(type=="clonal_normal") %>% group_by(project,patient) %>% tally(), by=c("project","patient"), all.x=TRUE) %>% 
  as_tibble() %>% 
  rename("patient_count"=n) %>% 
  group_by(patient) %>%
  add_count(VAR, name="VAR_patient_count") %>% 
  filter(VAR_patient_count > 1) %>% 
  select(-sample) %>% 
  group_by(patient,VAR) %>% 
  unique() %>% 
  ungroup()

  # edit variant count shared within the same late-branching lineages


 ## 2. calculate binomial probability ------------------------------------

df_shared_var_res <- df_shared_var %>% 
  filter(project=="Line1" | project=="DB" | project=="Hblood" | project=="blood_new" | project=="Abortus") %>% 
  filter(VAR!="414 T>G") %>% 
  arrange(patient,VAR) %>% 
  mutate(mle=(VAR_total_count-VAR_patient_count)/(total_count-patient_count)) %>% 
  rowwise() %>% 
  mutate(prob=sum(dbinom(x=VAR_patient_count_mod:patient_count_mod, size=patient_count_mod, prob=mle))) %>% 
  ungroup() %>% 
  mutate(vartype=ifelse(prob>=0.01, "recurrent", "fe"))