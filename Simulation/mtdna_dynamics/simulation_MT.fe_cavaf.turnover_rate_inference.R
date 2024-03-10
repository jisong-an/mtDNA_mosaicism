################################
#
#   HetFE turnover inference
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
#################################

## Inference of mitotic & homeostatic turnover rate



## Mitotic turnover ========================================================

load("~/R/R_script/11_Clone_MT/simul_fixvaf_result_fevar_sample_231229.RData")
library(parallel)

simul_result_fixvaf <- setClass("simul_result_fixvaf",
                                slots = c(
                                  pseudoVAF = "numeric",
                                  VAR = "character",
                                  patient = "character",
                                  tissue = "character",
                                  age = "numeric",
                                  fixvaf_simul_original = "data.frame",
                                  fixvaf_simul_original_dist = "data.frame",
                                  fixvaf_simul_sample = "data.frame",
                                  fixvaf_simul_sample_dist = "data.frame",
                                  obs = "data.frame"
                                ))

fixvaf_merge_sample <- do.call(rbind, mclapply(simul_result_fixvaf_list, function(x) {
  x@fixvaf_simul_sample_dist <- x@fixvaf_simul_sample_dist %>% slice(1:100) %>% 
    mutate(patient=x@patient, VAR=x@VAR, age=x@age, tissue=x@tissue, pseudoVAF=x@pseudoVAF)
  return(x@fixvaf_simul_sample_dist)
}, mc.cores=2))

fixvaf_merge_sample_mod <- fixvaf_merge_sample %>% 
  mutate(across(c(generation:distance,age,pseudoVAF), as.numeric)) %>% 
  mutate(age=ifelse(age==0,0.37,age)) %>%  # 19weeks = 0.37 yr
  mutate(gen_yr=generation/age) %>% 
  select(patient,VAR,age,tissue,pseudoVAF,gen_yr, distance)

fixvaf_merge_sample_mod2 <- fixvaf_merge_sample_mod %>% 
  mutate(group=paste0(patient,":",VAR)) %>% group_by(group) %>% 
  arrange(distance) %>% 
  slice(1:50) %>% 
  mutate(max=max(gen_yr), min=min(gen_yr),med=mean(gen_yr)) %>% 
  select(patient,VAR,age,tissue,pseudoVAF,group,max,min,med) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(genyr=paste(seq(round(min),round(max)),collapse=":")) %>% 
  ungroup() %>% 
  separate_rows(genyr,sep=":") %>% 
  select(patient,VAR,pseudoVAF,tissue,age,med,genyr) %>% 
  mutate(across(genyr,as.numeric)) %>% 
  filter(age>1) %>% 
  filter(genyr < 100) 


# peak
# colon
which.max(density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="colon"),]$genyr)$y) #166
density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="colon"),]$genyr)$x[166]  #14.28

# blood
which.max(density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="blood"),]$genyr)$y) #176
density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="blood"),]$genyr)$x[176]  #17.86

# fibroblast
which.max(density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="fibroblast"),]$genyr)$y) #150
density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="fibroblast"),]$genyr)$x[150]  #20.77


# graph: infer turnover
fixvaf_merge_sample_mod %>%
  mutate(group=paste0(patient,":",VAR)) %>% group_by(group) %>% 
  arrange(distance) %>% 
  slice(1:50) %>% 
  mutate(max=max(gen_yr), min=min(gen_yr),med=mean(gen_yr)) %>%
  mutate(group=paste0(age,":",patient,":",tissue)) %>%
  mutate(group2=paste0(age,":",patient,":",VAR)) %>%
  select(-gen_yr,-distance) %>% unique() %>% 
  filter(age>1) %>% 
  ggplot() +
  geom_errorbar(aes(x=reorder(group2,age),y=med, ymin=min, ymax=max), color="gray50", width=0.25) +
  geom_point(aes(x=reorder(group2,age),y=med, col=patient),size=2)  +
  facet_grid(~factor(tissue,levels=c("colon","fibroblast","blood")), scales="free_x",space="free") +
  theme_light()+
  xlab("fertilized-egg origin variants") +
  ylab("turnover per year") +
  theme(strip.text.x = element_text(color="black"), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(ylim=c(0,75)) + geom_hline(yintercept=c(14.28,17.86,20.77),lty=2,color="gray80")



## 95%CI interval - bootstrap
nBoot=1000
for (x in c("colon","fibroblast","blood")){
  
  mean.boot <- rep(NA,nBoot)
  tmp <- fixvaf_merge_sample_mod2 %>% filter(tissue==x) %>% arrange(genyr)  %>% pull(genyr)
  N <- length(tmp)
  
  for (i in 1:nBoot){
    temp_list <- sample(tmp,size=N,replace=TRUE)
    yval <- which.max(density(temp_list)$y)
    mean.boot[i] <- density(temp_list)$x[yval]
  }
  
  assign(paste0("mean.boot.",x),mean.boot)
  
}


# colon
quantile(mean.boot.colon,probs=c(0.025,0.975))
#  2.5%    97.5% 
#  9.72  17.55
quantile(mean.boot.fibroblast,probs=c(0.025,0.975))
#  2.5%    97.5% 
#  17.72  24.35
quantile(mean.boot.blood,probs=c(0.025,0.975))
#  2.5%    97.5% 
#  13.63  22.14

save(simul_result_fixvaf_list,fixvaf_merge_sample,fixvaf_merge_sample_mod,fixvaf_merge_sample_mod2,mean.boot.colon,mean.boot.fibroblast,mean.boot.blood, file="~/R/R_script/11_Clone_MT/simul_fe_fixvaf_mle_result.240125.RData" )





## Homeostatic turnover ====================================================

load("~/R/R_script/11_Clone_MT/simul_fixvaf_result_fevar_sample_231221.RData")
library(tidyverse)
library(parallel)

simul_result_fixvaf <- setClass("simul_result_fixvaf",
                                slots = c(
                                  pseudoVAF = "numeric",
                                  VAR = "character",
                                  patient = "character",
                                  tissue = "character",
                                  age = "numeric",
                                  fixvaf_simul_original = "data.frame",
                                  fixvaf_simul_original_dist = "data.frame",
                                  fixvaf_simul_sample = "data.frame",
                                  fixvaf_simul_sample_dist = "data.frame",
                                  obs = "data.frame"
                                ))


# merge simulation results of all samples
fixvaf_merge_sample <- do.call(rbind, mclapply(simul_result_fixvaf_list, function(x) {
  x@fixvaf_simul_sample_dist <- x@fixvaf_simul_sample_dist %>% slice(1:100) %>% 
    mutate(patient=x@patient, VAR=x@VAR, age=x@age, tissue=x@tissue, pseudoVAF=x@pseudoVAF)
  return(x@fixvaf_simul_sample_dist)
}, mc.cores=2))

fixvaf_merge_sample_mod <- fixvaf_merge_sample %>% 
  mutate(across(c(generation:distance,age,pseudoVAF), as.numeric)) %>% 
  mutate(age=ifelse(age==0,0.37,age)) %>%  # 19weeks = 0.37 yr
  mutate(gen_yr=generation/age) %>% 
  select(patient,VAR,age,tissue,pseudoVAF,gen_yr, distance)

fixvaf_merge_sample_mod2 <- fixvaf_merge_sample_mod %>% 
  mutate(group=paste0(patient,":",VAR)) %>% 
  group_by(group) %>% 
  arrange(distance) %>% 
  slice(1:30) %>% 
  mutate(max=max(gen_yr), min=min(gen_yr),med=mean(gen_yr)) %>% 
  select(patient,VAR,age,tissue,pseudoVAF,group,max,min,med) %>% 
  unique() %>% 
  rowwise() %>% 
  mutate(genyr=paste(seq(round(min),round(max)),collapse=":")) %>% 
  ungroup() %>% 
  separate_rows(genyr,sep=":") %>% 
  select(patient,VAR,pseudoVAF,tissue,age,med,genyr) %>% 
  mutate(across(genyr,as.numeric)) %>% 
  filter(age>1) %>% 
  filter(genyr < 100) 


# peak
# colon
which.max(density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="colon"),]$genyr)$y) #178
density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="colon"),]$genyr)$x[178]  #6.48

# blood
which.max(density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="blood"),]$genyr)$y) #198
density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="blood"),]$genyr)$x[198]  #9.44

# fibroblast
which.max(density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="fibroblast"),]$genyr)$y) #138
density(fixvaf_merge_sample_mod2[which(fixvaf_merge_sample_mod2$tissue=="fibroblast"),]$genyr)$x[138]  #11.53


# graph: infer turnover
fixvaf_merge_sample_mod %>%
  mutate(group=paste0(patient,":",VAR)) %>% group_by(group) %>% 
  arrange(distance) %>% 
  slice(1:50) %>% 
  mutate(max=max(gen_yr), min=min(gen_yr),med=mean(gen_yr)) %>%
  mutate(group=paste0(age,":",patient,":",tissue)) %>%
  mutate(group2=paste0(age,":",patient,":",VAR)) %>%
  select(-gen_yr,-distance) %>% unique() %>% 
  filter(age>1) %>% 
  ggplot() +
  geom_errorbar(aes(x=reorder(group2,age),y=med, ymin=min, ymax=max), color="gray50", width=0.25) +
  geom_point(aes(x=reorder(group2,age),y=med, col=patient),size=2)  +
  facet_grid(~factor(tissue,levels=c("colon","fibroblast","blood")), scales="free_x",space="free") +
  theme_light()+
  xlab("fertilized-egg origin variants") +
  ylab("turnover per year") +
  theme(strip.text.x = element_text(color="black"), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  coord_cartesian(ylim=c(0,75)) + geom_hline(yintercept=c(6.48,9.44,11.53),lty=2,color="gray80")


## 95%CI interval - bootstrap
nBoot=1000
for (x in c("colon","fibroblast","blood")){
  
  mean.boot <- rep(NA,nBoot)
  tmp <- fixvaf_merge_sample_mod2 %>% filter(tissue==x) %>% arrange(genyr)  %>% pull(genyr)
  N <- length(tmp)
  
  for (i in 1:nBoot){
    temp_list <- sample(tmp,size=N,replace=TRUE)
    yval <- which.max(density(temp_list)$y)
    mean.boot[i] <- density(temp_list)$x[yval]
  }
  
  assign(paste0("mean.boot.",x),mean.boot)
  
}


# colon
quantile(mean.boot.colon,probs=c(0.025,0.975))
#  2.5%    97.5% 
#  3.84  9.50
quantile(mean.boot.fibroblast,probs=c(0.025,0.975))
#  2.5%    97.5% 
#  9.68 14.23
quantile(mean.boot.blood,probs=c(0.025,0.975))
#  2.5%    97.5% 
#  6.35 12.35


save(simul_result_fixvaf_list,fixvaf_merge_sample,fixvaf_merge_sample_mod,fixvaf_merge_sample_mod2,mean.boot.colon,mean.boot.fibroblast,mean.boot.blood, file="~/R/R_script/11_Clone_MT/simul_fe_fixvaf_mle_result.mod1.240126.RData" )
