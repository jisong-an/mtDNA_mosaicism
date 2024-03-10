################################
#
#   mutation rate inference
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################


library(tidyverse)
load("~/R/R_script/11_Clone_MT/simul_result_mrate_mod2_240129.RData")
simul_result_mrate <- setClass("simul_result_mrate",
                               slots = c(
                                 patient = "character",
                                 tissue = "character",
                                 age = "numeric",
                                 mrate_simul = "data.frame",
                                 mrate_simul_dist = "data.frame",
                                 obs = "data.frame"
                               ))

simul_result_mrate_lst_mod2 <- simul_result_mrate_lst
mrate_dist_merge_mod2 <- mrate_dist_merge


## MLE ------------------------------------

calculate_distance_vaf_v2 <- function(simul,obs){
  
  simul_res <- simul %>%
    mutate(distance=((`5`-obs[["5"]])^2+
                       (`10`-obs[["10"]])^2+
                       (`15`-obs[["15"]])^2+
                       (`20`-obs[["20"]])^2+
                       (`25`-obs[["25"]])^2+
                       (`30`-obs[["30"]])^2+
                       (`35`-obs[["35"]])^2+
                       (`40`-obs[["40"]])^2+
                       (`45`-obs[["45"]])^2+
                       (`50`-obs[["50"]])^2+
                       (`55`-obs[["55"]])^2+
                       (`60`-obs[["60"]])^2+
                       (`65`-obs[["65"]])^2+
                       (`70`-obs[["70"]])^2+
                       (`75`-obs[["75"]])^2+
                       (`80`-obs[["80"]])^2+
                       (`85`-obs[["85"]])^2+
                       (`90`-obs[["90"]])^2+
                       (`95`-obs[["95"]])^2+
                       (`100`-obs[["100"]])^2+
                       (over3N-obs[["over3N"]])^2+
                       #(var0N-obs[["var0N"]])^2+
                       (over2N-obs[["over2N"]])^2))
  
  invisible(simul_res)
  
}

## change the method of distance calculation
mrate_dist_merge_mod2_v2 <- tibble()
for (i in 1:31){
  dist_temp <- calculate_distance_vaf_v2(simul_result_mrate_lst_mod2[[i]]@mrate_simul, simul_result_mrate_lst_mod2[[i]]@obs)
  dist_temp <- dist_temp %>% mutate(patient=simul_result_mrate_lst_mod2[[i]]@patient,tissue=simul_result_mrate_lst_mod2[[i]]@tissue, age=simul_result_mrate_lst_mod2[[i]]@age)
  mrate_dist_merge_mod2_v2 <- rbind(mrate_dist_merge_mod2_v2,dist_temp)
}

mrate_dist_merge_mod2_v2_uniq <- mrate_dist_merge_mod2_v2 %>% 
  group_by(patient) %>% 
  arrange(distance) %>% 
  dplyr::slice(1:50)%>% 
  mutate(mrate_min=min(mrate),mrate_max=max(mrate), mrate_avg=mean(mrate)) %>% 
  select(patient,age,tissue,mrate_min,mrate_max,mrate_avg) %>%
  unique() 

# graph
mrate_dist_merge_mod2_v2_uniq %>% 
  mutate(across(age,~replace(.,is.na(.),0)), across(tissue,~replace(.,is.na(.),"colon"))) %>% 
  arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% ungroup() %>% 
  mutate(num=row_number()) %>% 
  ggplot(aes(x=reorder(patient,num),y=mrate_avg)) + 
  geom_errorbar(aes(ymin=mrate_min,ymax=mrate_max),col="gray50",width=0.15) +
  geom_point(aes(col=tissue),size=2)  + 
  coord_cartesian(ylim=c(0,4e-7)) + 
  theme_classic() + 
  scale_color_manual(values = pal_mt) + xlab("patient") + ylab("mutation rate / bp / generation")


# peak
mrate_dist_merge_mod2_v2_uniq_forpeak <- mrate_dist_merge_mod2_v2_uniq %>% 
  rowwise() %>% 
  mutate(mrate_min=mrate_min*10^9,mrate_max=mrate_max*10^9) %>% 
  mutate(mrate=paste(seq(round(mrate_min),round(mrate_max)),collapse=":")) %>% 
  ungroup() %>% 
  separate_rows(mrate,sep=":") %>% 
  select(patient,age,tissue,mrate_avg,mrate) %>% 
  mutate(across(mrate,as.numeric)) %>% 
  filter(age>1) %>% 
  mutate(mrate=mrate/(10^9))


# peak
which.max(density(mrate_dist_merge_mod2_v2_uniq_forpeak$mrate)$y) #72
density(mrate_dist_merge_mod2_v2_uniq_forpeak$mrate)$x[72]  # 5e-8


save(simul_result_mrate,simul_result_mrate_lst_mod2,mrate_dist_merge_mod2,mrate_dist_merge_mod2_v2,mrate_dist_merge_mod2_v2_uniq,mrate_dist_merge_mod2_v2_uniq_forpeak, file="~/R/R_script/11_Clone_MT/simul_result_mrate_mod2_parse_240130.RData")


