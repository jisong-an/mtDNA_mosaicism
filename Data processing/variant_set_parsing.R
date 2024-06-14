################################
#
#   Variant set parsing
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

library(tidyverse)
 
getAbsolutePath <- function(dir, pattern) {    
  return(paste0(normalizePath(dir), "/", list.files(normalizePath(dir), pattern)))  
}


df <- tibble(project=character(), patient=character(), sample=character(), VAR=character(), POS=double(), REF  =character(), ALT=character(), Consequence=character(), Gene=character(), VAF=double())


# 1.normal colon ------------------------------------------------------------------
  
    file_line1 <- getAbsolutePath("/home/users/anjisong/project/11_Clone_MT/07_VarMat/01_Line1/Estcutoff/result_tier1_tier2/modify","*wogerm.tsv")
  
    for (f in file_line1) {
      df_temp <- read_tsv(f) %>% rename("patient"=sample) %>% 
        gather(16:(ncol(.)-6),key="sample", value="VAF") %>% 
        select(project,patient, sample,VAR,POS,REF,ALT,Consequence,Gene,VAF) %>% 
        filter(VAF>0)
      df <- rbind(df,df_temp)
      
    }
  
  
# 2. fibroblast ---------------------------------------------------------------------
  
    file_db <- getAbsolutePath("/home/users/anjisong/project/11_Clone_MT/07_VarMat/02_DB/Estcutoff/result_tier1_tier2/modify","*wogerm.tsv")
  
    for (f in file_db) {
      df_temp <- read_tsv(f) %>% rename("patient"=sample) %>% 
        gather(16:(ncol(.)-6),key="sample", value="VAF") %>% 
        select(project,patient, sample,VAR,POS,REF,ALT,Consequence,Gene,VAF) %>% 
        filter(VAF>0)
      df <- rbind(df,df_temp)
    }


# 3. HSC (1) ---------------------------------------------------------------------
  
    file_hblood <- getAbsolutePath("/home/users/anjisong/project/11_Clone_MT/07_VarMat/04_H.Blood/Estcutoff/result_tier1_tier2/modify","*wogerm.tsv")
  
    for (f in file_hblood) {
      df_temp <- read_tsv(f) %>% rename("patient"=sample) %>% 
        gather(16:(ncol(.)-6),key="sample", value="VAF") %>% 
        select(project,patient, sample,VAR,POS,REF,ALT,Consequence,Gene,VAF) %>% 
        filter(VAF>0)
      df <- rbind(df,df_temp)
    }


# 4. HSC (2) ---------------------------------------------------------------------
    
    file_blood_new <- getAbsolutePath("/home/users/anjisong/project/11_Clone_MT/07_VarMat/09_blood_new/Estcutoff/result_tier1_tier2/modify","*wogerm.tsv")
    
    for (f in file_blood_new) {
      df_temp <- read_tsv(f) %>% rename("patient"=sample) %>% 
        gather(16:(ncol(.)-6),key="sample", value="VAF") %>% 
        select(project,patient, sample,VAR,POS,REF,ALT,Consequence,Gene,VAF) %>% 
        filter(VAF>0)
      df <- rbind(df,df_temp)
    }
  

# 5. normal colon (abortus) -------------------------------------------------------------------------
    
    df_temp <- read_tsv("~/project/11_Clone_MT/07_VarMat/11_Abortus/Estcutoff/result_tier1_tier2/modify/Abortus_SA03_VarMat.filt.nb0.rescue.tier2.tier1_2.wogerm.tsv")%>% 
      rename("patient"=sample) %>% 
      gather(16:(ncol(.)-6),key="sample", value="VAF") %>% 
      select(project,patient, sample,VAR,POS,REF,ALT,Consequence,Gene,VAF) %>% 
      filter(VAF>0)
    df <- df %>% rbind(df_temp)
  

# 6. MUTYH adenoma -------------------------------------------------------------------------
    
    df_temp <- read_tsv("~/project/11_Clone_MT/07_VarMat/10_MUTYH_colon/Estcutoff/result_tier1_tier2/modify/MUTYH_HC22_VarMat.filt.nb0.rescue.tier2.tier1_2.wogerm.tsv")%>% 
      rename("patient"=sample) %>% 
      gather(16:(ncol(.)-6),key="sample", value="VAF") %>% 
      select(project,patient, sample,VAR,POS,REF,ALT,Consequence,Gene,VAF) %>% 
      filter(VAF>0)
    df <- df %>% rbind(df_temp)  
    
    
# 7. Merge with Sample info ------------------------------------------------
  
    df <- df %>% 
      merge(sample_info_final_v3 %>% select(1:7), by=c("project","patient","sample"), all=TRUE) %>% 
      as_tibble() 
