################################
#
#   Figure 3 -Lifetime drift
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

library(tidyverse)
library(ggsci)
library(ggtree)
library(RColorBrewer)
library(grid)
library(circlize)


### Fig3a. HetFE VAF distribution -----------------------------

# all variant (with pseudoVAF)
fe_var_df <- df %>% filter(vartype=="fe") %>% 
  filter(type=="clonal_normal") %>%
  group_by(patient, VAR) %>% 
  mutate(min_vaf=min(VAF),max_vaf=max(VAF)) %>% 
  ungroup() %>% 
  arrange(factor(tissue2, levels=c("colon","fibroblast","blood")), age, desc(pseudoVAF)) %>%
  mutate(numbering=row_number()) %>% 
  select(patient,VAR,sample,VAF,age,tissue2, patientN, pseudoVAF, numbering, min_vaf, max_vaf) %>% 
  mutate(group=paste0(patient,":",VAR)) 

# graph
pat_pos <- fe_var_df %>% select(patient,group) %>% unique() %>% group_by(patient) %>% mutate(n=n()) %>% ungroup() %>% select(patient,n) %>% unique() %>% pull(n)
pat_age <- fe_var_df %>% select(patient,patientN, age) %>% unique() %>% pull(age)
pat_clone <- fe_var_df %>% select(patient,patientN, age) %>% unique() %>% pull(patientN)
pat_tissue <- fe_var_df %>% select(patient,patientN, age, tissue2) %>% unique() %>% pull(tissue2)

# normalized pseudo VAF
pat_pvaf <- fe_var_df %>% select(patient,patientN, age, tissue2,pseudoVAF) %>% unique() %>% 
  mutate(logp=log2(pseudoVAF)) %>% 
  mutate(up=ifelse(2**floor(logp)>0.125,2**ceiling(logp),0.125), down=ifelse(2**floor(logp)>0.125,2**floor(logp),0)) %>% 
  mutate(intv=up-down) %>% 
  mutate(pseudoVAF_norm=(pseudoVAF-down)/intv) %>% 
  mutate(pseudoVAF_norm_final = ifelse(2**floor(logp)>0.125,floor(logp)+pseudoVAF_norm,pseudoVAF_norm)) %>% 
  mutate(pseudoVAF_norm_final_plot = ifelse(2**floor(logp)>0.125,pseudoVAF_norm_final+3,pseudoVAF_norm_final))%>% 
  pull(pseudoVAF_norm_final_plot)


xaxis_pos <- c(0)
for (i in 1:length(pat_pos)){xaxis_pos[i+1]=sum(pat_pos[1:i])}
xaxis_plot_pos = xaxis_pos/sum(pat_pos)

grid.newpage()

pushViewport(viewport(x=0,y=0, width=1,height=1,just=c('left','bottom')))


## tissue [0.9-0.95]
pushViewport(viewport(x=0,y=0.9, height=0.05, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

grid.text('tissue', x=-0.03, gp=gpar(fontsize=8))

for (i in 1:length(pat_clone)){
  tissue = pat_tissue[i]
  if (tissue=="colon") {bgcolor=pal_mt[1]}
  else if (tissue=="fibroblast") { bgcolor=pal_mt[2]}
  else{bgcolor=pal_mt[3]}
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
  #grid.text(pat_name[i], x=(xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=-0.04, gp=gpar(fontsize=7.5))
}
popViewport(3)


# clone count [0.85-0.9]
pushViewport(viewport(x=0,y=0.85, height=0.05, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

grid.text('# of clone', x=-0.03, gp=gpar(fontsize=8))
col_fun = colorRamp2(c(10,50,400),c("white",pal_jco()(10)[2],pal_jama()(10)[2]))
#col_fun = colorRamp2(c(0,50,400),c(pal_jco()(10)[2],"white",pal_jama()(10)[2]))
#col_fun = colorRamp2(c(0,30,100),c("#FDFBF5","#FF97B5","#6E2A85"))

for (i in 1:length(pat_pos)){
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=col_fun(pat_clone[i]), col=FALSE))
  grid.text(pat_clone[i], x=(xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=0.5, gp=gpar(fontsize=6.5))
}
popViewport(2)



# VAF dist [0.2-0.85]
pushViewport(viewport(x=0, y=0.2, width=1, height=0.65, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

# background
for (i in 1:length(pat_pos)){
  if(i %% 2 == 0){bgcolor="gray92"}
  else{bgcolor = FALSE}
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
}

grid.yaxis(at=c(0,0.25,0.5,0.75,1), label=c(0,0.25,0.5,0.75,1), gp=gpar(fontsize=7))
grid.xaxis(label = FALSE, at=c(0,1))

grid.text('VAF', x=-0.03, rot=90, gp=gpar(fontsize=9))
grid.text('fertilized-egg origin variants', y=-0.04, gp=gpar(fontsize=9))

grid.segments(x0=xaxis_plot_pos[c(16,22)], x1=xaxis_plot_pos[c(16,22)],y0=0,y1=1, gp=gpar(color="gray60", fill=FALSE, lwd=0.5))

gt <- fe_var_df %>% ggplot(aes(x=reorder(group,numbering),y=VAF)) +
  geom_errorbar(aes(ymin=0,ymax=max_vaf), width=0.01, color="gray50") + 
  #geom_point(aes(y=pseudoVAF*2), col="gray30", fill="gray30", size=1.5, shape=23) +
  geom_point(aes(col=tissue2)) + 
  #scale_y_continuous(name="VAF (%)", sec.axis=sec_axis(~.*0.5, name="pseudoVAF (%)")) + 
  scale_color_manual(values = pal_mt) + 
  theme_void() + coord_cartesian(ylim=c(0,97))+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  xlab("fertilized-egg origin variant (order : pseudoVAF)")
g <- ggplotGrob(gt)
grid.draw(g)

popViewport(2)


## pseudoVAF [0.05-0.2]
pushViewport(viewport(x=0,y=0.05, height=0.15, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

grid.yaxis(at=seq(0,1,0.1), label=c(64,32,16,8,4,2,1,0.5,0.25,0.125,0), gp=gpar(fontsize=7))
grid.xaxis(label = FALSE, at=c(0,1))
grid.text('HEF', x=-0.03, rot=90, gp=gpar(fontsize=9))


# background
for (i in 1:length(pat_pos)){
  if(i %% 2 == 0){bgcolor="gray92"}
  else{bgcolor = FALSE}
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
}

# y-value : max VAF
xaxis_pos <- seq(1,length(pat_pvaf))/length(pat_pvaf)
y_value_vaf <- pat_pvaf/9

# bar graph
for (i in 1:length(pat_pvaf)){
  grid.rect(x=xaxis_pos[i],y=1, just=c('right','top'), height=y_value_vaf[i], width=xaxis_data[1], gp=gpar(fill="gray30", col="gray", lwd=0.25))
}
popViewport(2)

## additional age
# graph
pat_pos <- fe_var_df %>% select(patient,group) %>% unique() %>% group_by(patient) %>% mutate(n=n()) %>% ungroup() %>% select(patient,n) %>% unique() %>% pull(n)
pat_age <- fe_var_df %>% select(patient,patientN, age) %>% unique() %>% pull(age)

xaxis_pos <- c(0)
for (i in 1:length(pat_pos)){xaxis_pos[i+1]=sum(pat_pos[1:i])}
xaxis_plot_pos = xaxis_pos/sum(pat_pos)

grid.newpage()

pushViewport(viewport(x=0,y=0, width=1,height=1,just=c('left','bottom')))

## age [0.9-0.95]
pushViewport(viewport(x=0,y=0.9, height=0.05, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

grid.text('age', x=-0.03, gp=gpar(fontsize=8))
col_fun = colorRamp2(c(0,30,100),c("lightgoldenrod","goldenrod1","darkorchid4"))
#col_fun = colorRamp2(c(0,30,100),c("#FDFBF5","#FF97B5","#6E2A85"))

for (i in 1:length(pat_pos)){
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=col_fun(pat_age[i]), col=FALSE))
  #grid.text(pat_name[i], x=(xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=-0.04, gp=gpar(fontsize=7.5))
}
popViewport(3)



### Fig3b. ideal VAF & real VAF histogram (16256 C>T) -----------------------

# DB10 16256 C>T
df %>% filter(patient=="10" & POS==16256) %>% 
  select(patient,sample,VAF) %>% 
  merge(sample_info_final_v3 %>% 
          filter(patient=="10") %>% 
          select(patient,sample),by=c("patient","sample"), all=TRUE) %>% 
  as_tibble()  %>% 
  mutate(across(VAF,~replace(.,is.na(.),0))) %>% 
  mutate(VAF=VAF/100) %>% 
  ggplot() + 
  geom_histogram(aes(x=VAF), fill=pal_mt[2], color="white")  + 
  geom_vline(xintercept=0.32, color="gray30", lty=2) + 
  geom_hline(yintercept=0, color="gray90") +  
  theme_classic()



### Fig3d. abortus vaf & older --------------------------------

library(ggridges)

fe_vaf_all %>% filter(VAR %in% c("15649 A>G","12359 C>T","6378 T>C","8416 C>T"))%>% 
  ggplot(aes(x=VAF/100,y=VAR, fill=VAR)) + 
  geom_density_ridges( scale=1.25,alpha=0.8) + 
  theme_ridges() + 
  scale_fill_jco() + 
  theme(axis.title.y = element_blank(), axis.title.x = element_text(hjust = 0.5), legend.title = element_blank())  + 
  xlab("VAF")

# A tibble: 4 × 4
# patient VAR       pseudoVAF   age
# 3       15649 A>G      5.99    49
# 9       12359 C>T      6.09    82
# HC08    6378 T>C       7.18    53
# SA03    8416 C>T       6.33     0



### Fig3e. age-fixation index ----------------------------

fe_var_df_uniq <- fe_var_df %>% 
  filter(type=="clonal_normal")   %>% select(patient,VAR,age,tissue2, pseudoVAF) %>% unique()

fe_temp <- tibble()

for (i in 1:nrow(fe_var_df_uniq)){
  
  pat <- fe_var_df_uniq[[i,"patient"]]
  var <- fe_var_df_uniq[[i,"VAR"]]
  
  temp <- fe_var_df %>% 
    filter(type=="clonal_normal") %>% 
    filter(patient==pat) %>% 
    filter(VAR==var) %>% 
    select(patient,sample,VAR, VAF) %>% 
    merge(sample_info_final_v3 %>% 
            filter(patient==pat) %>% filter(type=="clonal_normal") %>%  
            select(patient,sample), by=c("patient","sample"), all=TRUE) %>% 
    as_tibble() %>% 
    mutate(across(VAR, ~replace(., is.na(.), var)), across(VAF, ~replace(., is.na(.),0))) 
  fe_temp <- rbind(fe_temp, temp)
}

fe_vaf_all <- fe_temp %>%  merge(fe_var_df_uniq, by=c("patient","VAR")) %>% as_tibble()


# fixation index
fe_vaf_all %>% 
  mutate(VAF=VAF/100, pseudoVAF=pseudoVAF/100) %>% 
  group_by(patient,VAR) %>% 
  mutate(VAFvar=var(VAF)) %>% 
  mutate(n=sum(VAF>0)) %>% 
  filter(n>1) %>% 
  select(patient,VAR,age,tissue2,pseudoVAF,VAFvar) %>% 
  unique() %>% 
  mutate(exp.freq=pseudoVAF*(1-pseudoVAF)) %>% 
  mutate(f_st=VAFvar/exp.freq) %>% 
  ungroup() %>% filter(pseudoVAF>0.01)  %>% 
  ggplot(aes(x=age, y=f_st)) + 
  geom_smooth(method='lm', level=0.95, fullrange=TRUE, color="gray50", fill="gray93") + 
  geom_point(aes(col=pseudoVAF), size=2) + 
  ylab("fixation index") + 
  theme_classic() + 
  scale_color_viridis(rescaler = function(x, to = c(0, 0.95), from = NULL) {
    ifelse(x<0.1, scales::rescale(x,to = to, from = c(min(x, na.rm = TRUE), 0.1)),0.95)}, option="inferno", direction=-1)

# cor.test
temp <- fe_vaf_all %>% 
  mutate(VAF=VAF/100, pseudoVAF=pseudoVAF/100) %>% 
  group_by(patient,VAR) %>% 
  mutate(VAFvar=var(VAF)) %>% 
  mutate(n=sum(VAF>0)) %>% 
  filter(n>1) %>% 
  select(patient,VAR,age,tissue2,pseudoVAF,VAFvar) %>% 
  unique() %>% 
  mutate(exp.freq=pseudoVAF*(1-pseudoVAF)) %>% 
  mutate(f_st=VAFvar/exp.freq) %>% 
  ungroup()


cor.test(temp %>% filter(pseudoVAF>0.01) %>% pull(age),temp %>% filter(pseudoVAF>0.01) %>% pull(f_st)) # pval=2.5*e-4, cor=0.640




### Fig3f. branching time --------------------

# 8th division = 3.8+1.2*7 ~= 12
# 16th division = 3.8+1.2*15 ~= 22
# 32th division = 3.8+1.2*31 ~= 41

load("~/R/R_script/11_Clone_MT/lineage_sharedvar.RData")

df_tree_shared_length_vaf_uniq %>% 
  filter(shared_length<42) %>% 
  mutate(group=ifelse(shared_length == 0, "A",
                      ifelse(shared_length <= 12, "B",
                             ifelse(shared_length <= 22, "C", "D")))) %>% 
  ggplot(aes(x=group,y=VAFdiff)) + 
  geom_point(aes(col=group), position=position_jitter(seed=1, width=0.2, height=0.2)) + 
  geom_violin(scale="width", alpha=0, width=0.8) +
  theme_classic() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(10,"BuPu"))(20)[c(6,9,12,16)]) +
  xlab("Branching time of two clones") + 
  ylab("Difference in heteroplasmy level of FE mutations") + theme(legend.position = "none")




###Fig3g. turnover count to fixation (mitotic) -----------------------------------------------

fixvaf_fixgen_df_mod_240122 %>% 
  filter(model=="model2" & mtCN==750) %>% 
  filter(fixVAF==1) %>% 
  ggplot() + 
  geom_point(aes(x=initVAF,y=meanGen, col=as.character(fixVAF)),size=2) + 
  theme_classic() + 
  xlab("caVAF") + 
  ylab("Avg. turnover for homoplasmy") + 
  theme(axis.title = element_text(size=10), axis.text = element_text(size=8)) + 
  ylim(c(200,1600)) + 
  scale_color_manual(values=c(pal_cosmic("hallmarks_light")(10)[2])) + 
  theme(legend.position = "none")



### Fig3h. compare results of mitotic & homeostatic turnover model --------------------------

fixvaf_fixgen_df_mod_240122 %>% 
  filter(mtCN==750) %>% 
  filter(fixVAF==1) %>% 
  select(fixVAF,meanGen,initVAF,mtCN,model) %>% 
  spread(key="model",value=meanGen) %>% 
  ggplot(aes(x=model2,y=model1)) + 
  geom_point() + 
  geom_smooth()



### Fig3i. turnover rate inference in two models --------------------------------------------


# mitotic turnover

load("~/R/R_script/11_Clone_MT/simul_fe_fixvaf_mle_result.240125.RData")


res_mle_fixvaf_mod2 <- fixvaf_merge_sample_mod_240125 %>%
  mutate(group=paste0(patient,":",VAR)) %>% group_by(group) %>% 
  arrange(distance) %>% 
  slice(1:50) %>% 
  mutate(max=max(gen_yr), min=min(gen_yr),med=mean(gen_yr)) %>%
  mutate(group=paste0(age,":",patient,":",tissue)) %>%
  mutate(group2=paste0(age,":",patient,":",VAR)) %>%
  select(-gen_yr,-distance) %>% unique() %>% 
  filter(age>1)


# graph
pat_pos <- res_mle_fixvaf_mod2 %>% arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age)%>% select(patient) %>% group_by(patient) %>% mutate(n=n()) %>% ungroup() %>% unique() %>% pull(n)
xaxis_pos <- c(0)
for (i in 1:length(pat_pos)){xaxis_pos[i+1]=sum(pat_pos[1:i])}
xaxis_plot_pos = xaxis_pos/sum(pat_pos)

grid.newpage()
pushViewport(viewport(x=0.1, y=0.1, width=0.65, height=0.5, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

# background
for (i in 1:length(pat_pos)){
  if(i %% 2 == 0){bgcolor="gray92"}
  else{bgcolor = FALSE}
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
}


# peak & 95% CI 
c_peak = 14.28
c_low = 9.72
c_high = 17.55
f_peak = 20.77
f_low = 17.72
f_high = 24.35
b_peak = 17.86
b_low = 13.63
b_high = 22.14
tot_height = 80


grid.rect(x=0,y=c_low/tot_height, just=c('left','bottom'), height=(c_high-c_low)/tot_height, width=xaxis_plot_pos[9], gp=gpar(col=FALSE, fill="gray80"))
grid.rect(x=xaxis_plot_pos[9],y=f_low/tot_height, just=c('left','bottom'), height=(f_high-f_low)/tot_height, width=xaxis_plot_pos[15]-xaxis_plot_pos[9], gp=gpar(col=FALSE, fill="gray80"))
grid.rect(x=xaxis_plot_pos[15],y=b_low/tot_height, just=c('left','bottom'), height=(b_high-b_low)/tot_height, width=1-xaxis_plot_pos[15], gp=gpar(col=FALSE, fill="gray80"))

grid.segments(x0=0,x1=xaxis_plot_pos[9], y0=c_peak/tot_height, y1=c_peak/tot_height, gp=gpar(color="gray30", fill=FALSE, lwd=0.5))
grid.segments(x0=xaxis_plot_pos[9],x1=xaxis_plot_pos[15], y0=f_peak/tot_height, y1=f_peak/tot_height, gp=gpar(color="gray30", fill=FALSE, lwd=0.5))
grid.segments(x0=xaxis_plot_pos[15],x1=1, y0=b_peak/tot_height, y1=b_peak/tot_height, gp=gpar(color="gray30", fill=FALSE, lwd=0.5))


grid.rect(gp=gpar(lty="solid", fill=FALSE))
grid.yaxis(at=c(0,0.25,0.5,0.75,1), label=c(0,20,40,60,80), gp=gpar(fontsize=7))
grid.text('turnover per year', x=-0.06, rot=90, gp=gpar(fontsize=9))
grid.text('fertilized-egg origin variants', y=-0.04, gp=gpar(fontsize=9))

grid.segments(x0=0,x1=1, y0=c(0.125,0.25,0.375,0.5,0.625,0.75,0.875), y1=c(0.125,0.25,0.375,0.5,0.625,0.75,0.875), gp=gpar(color="gray85", fill=FALSE, lwd=0.2))
grid.segments(x0=xaxis_plot_pos[c(9,15)], x1=xaxis_plot_pos[c(9,15)],y0=0,y1=1, gp=gpar(color="gray60", fill=FALSE, lwd=0.5))




gt <- res_mle_fixvaf_mod2 %>% mutate(group=paste0(tissue,":",age)) %>% 
  ggplot() + 
  geom_errorbar(aes(x=factor(group2,levels=c(res_mle_fixvaf_mod1 %>% arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% pull(group2))),
                    y=med, ymin=min, ymax=max), color="gray50", width=0.25) +
  geom_point(aes(x=factor(group2,levels=c(res_mle_fixvaf_mod1 %>% arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% pull(group2))),
                 y=med, col=tissue))  + 
  #facet_grid(~tissue, scales="free_x",space="free") + 
  coord_cartesian(ylim=c(3,75)) + 
  theme_void() + 
  theme(legend.position="none", strip.text = element_blank(), panel.spacing.x = unit(-0.2,"mm")) + 
  scale_color_manual(values=pal_mt)
g <- ggplotGrob(gt)
grid.draw(g)

popViewport(2)


# homeostatic turnover

load("~/R/R_script/11_Clone_MT/simul_fe_fixvaf_mle_result.mod1.240126.RData")


res_mle_fixvaf_mod1 <- fixvaf_merge_sample_mod_240126 %>%
  mutate(group=paste0(patient,":",VAR)) %>% group_by(group) %>% 
  arrange(distance) %>% 
  slice(1:50) %>% 
  mutate(max=max(gen_yr), min=min(gen_yr),med=mean(gen_yr)) %>%
  mutate(group=paste0(age,":",patient,":",tissue)) %>%
  mutate(group2=paste0(age,":",patient,":",VAR)) %>%
  select(-gen_yr,-distance) %>% unique() %>% 
  filter(age>1)


# graph
pat_pos <- res_mle_fixvaf_mod1 %>% arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age)%>% select(patient) %>% group_by(patient) %>% mutate(n=n()) %>% ungroup() %>% unique() %>% pull(n)
xaxis_pos <- c(0)
for (i in 1:length(pat_pos)){xaxis_pos[i+1]=sum(pat_pos[1:i])}
xaxis_plot_pos = xaxis_pos/sum(pat_pos)

grid.newpage()
pushViewport(viewport(x=0.1, y=0.1, width=0.65, height=0.5, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

# background
for (i in 1:length(pat_pos)){
  if(i %% 2 == 0){bgcolor="gray92"}
  else{bgcolor = FALSE}
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
}


# peak & 95% CI 
c_peak = 6.48
c_low = 3.84
c_high = 9.5
f_peak = 11.53
f_low = 9.68
f_high = 14.23
b_peak = 9.44
b_low = 6.35
b_high = 12.35
tot_height = 80


grid.rect(x=0,y=c_low/tot_height, just=c('left','bottom'), height=(c_high-c_low)/tot_height, width=xaxis_plot_pos[9], gp=gpar(col=FALSE, fill="gray80"))
grid.rect(x=xaxis_plot_pos[9],y=f_low/tot_height, just=c('left','bottom'), height=(f_high-f_low)/tot_height, width=xaxis_plot_pos[15]-xaxis_plot_pos[9], gp=gpar(col=FALSE, fill="gray80"))
grid.rect(x=xaxis_plot_pos[15],y=b_low/tot_height, just=c('left','bottom'), height=(b_high-b_low)/tot_height, width=1-xaxis_plot_pos[15], gp=gpar(col=FALSE, fill="gray80"))

grid.segments(x0=0,x1=xaxis_plot_pos[9], y0=c_peak/tot_height, y1=c_peak/tot_height, gp=gpar(color="gray30", fill=FALSE, lwd=0.5))
grid.segments(x0=xaxis_plot_pos[9],x1=xaxis_plot_pos[15], y0=f_peak/tot_height, y1=f_peak/tot_height, gp=gpar(color="gray30", fill=FALSE, lwd=0.5))
grid.segments(x0=xaxis_plot_pos[15],x1=1, y0=b_peak/tot_height, y1=b_peak/tot_height, gp=gpar(color="gray30", fill=FALSE, lwd=0.5))


grid.rect(gp=gpar(lty="solid", fill=FALSE))
grid.yaxis(at=c(0,0.25,0.5,0.75,1), label=c(0,20,40,60,80), gp=gpar(fontsize=7))
grid.text('turnover per year', x=-0.06, rot=90, gp=gpar(fontsize=9))
grid.text('fertilized-egg origin variants', y=-0.04, gp=gpar(fontsize=9))

grid.segments(x0=0,x1=1, y0=c(0.125,0.25,0.375,0.5,0.625,0.75,0.875), y1=c(0.125,0.25,0.375,0.5,0.625,0.75,0.875), gp=gpar(color="gray85", fill=FALSE, lwd=0.2))
grid.segments(x0=xaxis_plot_pos[c(9,15)], x1=xaxis_plot_pos[c(9,15)],y0=0,y1=1, gp=gpar(color="gray60", fill=FALSE, lwd=0.5))




gt <- res_mle_fixvaf_mod1 %>% mutate(group=paste0(tissue,":",age)) %>% 
  ggplot() + 
  geom_errorbar(aes(x=factor(group2,levels=c(res_mle_fixvaf_mod1 %>% arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% pull(group2))),
                    y=med, ymin=min, ymax=max), color="gray50", width=0.25) +
  geom_point(aes(x=factor(group2,levels=c(res_mle_fixvaf_mod1 %>% arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% pull(group2))),
                 y=med, col=tissue))  + 
  #facet_grid(~tissue, scales="free_x",space="free") + 
  coord_cartesian(ylim=c(3,75)) + 
  theme_void() + 
  theme(legend.position="none", strip.text = element_blank(), panel.spacing.x = unit(-0.2,"mm")) + 
  scale_color_manual(values=pal_mt)
g <- ggplotGrob(gt)
grid.draw(g)

popViewport(2)



