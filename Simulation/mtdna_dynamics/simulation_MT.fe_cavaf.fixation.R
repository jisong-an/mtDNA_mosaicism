################################
#
#   Average turnover to reach homoplasmy (fixation)
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
#################################


## load simulation result ----------------------------------

fixvaf_lst_fixgen <- vector(mode="list", length=14)
names(fixvaf_lst_fixgen) <- round(c(0.001, 0.005, 0.01, 0.03, 0.05, seq(0.1,0.9,0.1)),4)

file <- getAbsolutePath("~/project/11_Clone_MT/10_Simulation/231229_fixvaf/", "*fixdata.iter1.tsv")

for (i in 1:length(file)){
  
  t <- read_tsv(file[i])
  fixvaf_lst_fixgen[[i]] <- t
  
}


## histogram & mean generation -------------------------------------------

fixvaf_fixgen_df <- tibble()

for (i in 1:length(fixvaf_lst_fixgen)){
  
  vafname = names(fixvaf_lst_fixgen)[i]
  temp <- fixvaf_lst_fixgen[[i]] %>% 
    group_by(fixVAF) %>% 
    summarise(meanGen=mean(fixGen), medianGen=median(fixGen), 
              minGen=min(fixGen), maxGen=max(fixGen), sdGen=round(sd(fixGen),2)) %>% 
    merge(fixvaf_lst_fixgen[[i]] %>% 
            group_by(fixVAF) %>% 
            tally() %>% 
            mutate(n=n/10000) %>% 
            dplyr::rename("cell_ratio"=n), by=c("fixVAF")) %>% 
    as_tibble() %>% 
    mutate(initVAF=vafname) %>% 
    mutate(unfix=10000-nrow(fixvaf_lst_fixgen[[i]]))

  fixvaf_fixgen_df <- rbind(fixvaf_fixgen_df,temp)
}

save(fixvaf_lst_fixgen, fixvaf_fixgen_df, freq_list,file="~/R/R_script/11_Clone_MT/simul_fixvaf_result_lowMem.fixationTurnover.240109.RData")
