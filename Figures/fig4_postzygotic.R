################################
#
#   Figure 4 -post-zygotic mutation
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


### Fig4a. # of mut VS. age ----------------------------------

normal_mut <- df %>% 
  filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(type=="clonal_normal") %>% 
  filter(vartype!="fe" & vartype!="fe_gray") %>% 
  filter(!is.na(VAR)) %>% 
  group_by(project,patient,age,sample,type, tissue2) %>% 
  mutate(sumVAF=sum(VAF), n=n(), maxVAF=max(VAF)) %>% 
  group_by(project,patient,age,sample,type, tissue2,sumVAF,maxVAF) %>% 
  tally() %>% 
  merge(sample_info_final_v3 %>% filter(type=="clonal_normal")  %>% 
          filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
          select(project,patient,sample,age,type, tissue2),
        by=c("project","patient","sample","age","type","tissue2"),all=TRUE) %>%
  as_tibble() %>% 
  mutate(across(n,~replace(.,is.na(.),0)), across(sumVAF,~replace(.,is.na(.),0)), across(maxVAF, ~replace(.,is.na(.),0))) %>% 
  merge(mutsig_df, by=c("patient","sample"), all.x=TRUE) %>% 
  as_tibble() %>% 
  filter(SBS7num < 1500 | is.na(SBS7num))


# count
normal_mut %>% 
  group_by(project,patient,age,type,tissue2) %>%
  filter(age!=0) %>% 
  summarise(mean_n=mean(n), sd_n=sd(n),n=n()) %>% 
  mutate(q1=mean_n-1.96*sd_n/sqrt(n), q3=mean_n+1.96*sd_n/sqrt(n)) %>% 
  ggplot(aes(x=age,y=mean_n))  + 
  geom_errorbar(aes(ymin=q1, ymax=q3),col="gray", width=1) +
  geom_point(aes(col=tissue2),size=1.5) + 
  geom_smooth(method='lm',level=0.95, fullrange=TRUE, color="gray30", fill="lightgray") + 
  theme_classic() + 
  scale_color_manual(values=pal_mt) + 
  ggtitle("average mutation count per patient")


# cor.test
cor.test(normal_mut %>% filter(age!=0) %>% 
           group_by(project,patient,age,type,tissue2) %>%
           summarise(mean_n=mean(n), sd_n=sd(n),n=n()) %>% pull(age), normal_mut %>% filter(age!=0) %>% 
           group_by(project,patient,age,type,tissue2) %>%
           summarise(mean_n=mean(n), sd_n=sd(n),n=n()) %>% pull(mean_n))
# cor = 0.282, pval= 0.131



### Fig4b. maxVAF distribution per age ------------------------------------

df %>% 
  filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(type=="clonal_normal") %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(vartype!="fe") %>% 
  filter(vartype!="fe_gray") %>% 
  filter(!is.na(VAR)) %>% 
  merge(sample_info_final_v3 %>% 
          filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
          filter(type=="clonal_normal") %>% select(patient,sample,age,tissue2), by=c("patient","sample","age","tissue2"), all=TRUE) %>% 
  as_tibble() %>% 
  mutate(across(VAF, ~replace(., is.na(.), 0))) %>% 
  group_by(patient,sample,age, tissue2) %>% 
  summarise(max=max(VAF, na.rm=TRUE))  %>% 
  merge(mutsig_df, by=c("patient","sample"), all.x=TRUE) %>% 
  as_tibble() %>% 
  mutate(SBS7num=SBS7*SBSnum/100) %>% 
  mutate(group=paste0(age," (",patient, ")")) %>% 
  filter(SBS7num < 1500 | is.na(SBS7num)) %>% 
  mutate(group_vaf=ifelse(max>=90, "VAF>90", ifelse(max>=75, "VAF>75", ifelse(max>=50, "VAF>50", ifelse(max>=25, "VAF>25", ifelse(max>0,"VAF<25", "VAF<=0")))))) %>%  
  ungroup() %>% 
  group_by(group, group_vaf, tissue2) %>% 
  tally() %>% 
  group_by(group,tissue2) %>% 
  mutate(total=sum(n)) %>% 
  mutate(ratio=n/total*100) %>% 
  mutate(pos=cumsum(ratio)-(0.5*ratio)) %>% 
  rename(tissue2="tissue") %>% 
  ggplot() + 
  geom_bar(aes(fill=factor(tissue,levels=c("colon","fibroblast","blood")), x=group, y=n, group=group_vaf, alpha=group_vaf),position="fill", stat="identity") + 
  geom_text(aes(x=group, y=1-pos/100,label=round(ratio))) + 
  theme_light() + ylab("") +  
  ggtitle("maximum VAF per patient (only use UV-mediated variant < 1500 clones)") + 
  facet_grid(~factor(tissue,levels=c("colon","fibroblast","blood")), scales = "free_x", space = "free") + 
  scale_fill_manual(values=c(brewer.pal("Blues",n=6)[6],brewer.pal("Reds",n=6)[6],brewer.pal("Greens",n=6)[6])) + 
  scale_alpha_manual(values=c(0.05,0.2,0.4,0.6,0.8,1)) + 
  theme(axis.text.x = element_text(angle=25), panel.grid = element_blank(), strip.text.x = element_text(color="black",face="bold"), legend.position="top") + labs(fill="tissue")



### Fig4c. DB8 & HC10 ----------------------------------

# DB8 (1-3 mutations only)
db8_vaf_comp <- df %>% filter(patient=="8") %>% 
  filter(type=="clonal_normal") %>% 
  filter(vartype!="fe") %>% 
  filter(!(VAR %in% c(freq_recur_var_proj %>% filter(project=="DB") %>% pull(VAR)))) %>% 
  mutate(across(VAF, ~replace(.,is.na(.),0))) %>% 
  select(patient,VAR,sample,VAF)%>% 
  group_by(sample) %>% 
  mutate(maxVAF=max(VAF)) %>% 
  filter(maxVAF+VAF < 105 | maxVAF==VAF) %>% 
  filter(!(maxVAF>95 & maxVAF!=VAF)) %>% 
  arrange(sample,desc(VAF)) %>% 
  dplyr::slice(1:5)  %>% 
  mutate(secondVAF=sort(VAF,decreasing=TRUE)[2]) %>% 
  mutate(thirdVAF=sort(VAF,decreasing=TRUE)[3]) %>% 
  mutate(maxsecond=maxVAF+secondVAF) %>% 
  mutate(maxthird=maxVAF+secondVAF+thirdVAF) %>% 
  mutate(group=ifelse(maxsecond>105 | is.na(maxsecond), "1mut", ifelse(maxthird>105 | is.na(maxthird), "2mut","3mut"))) %>%
  mutate(sumVAF=ifelse(group=="1mut",maxVAF, ifelse(group=="2mut",maxsecond, maxthird))) %>% 
  dplyr::slice(1:as.numeric(gsub("mut","",group))) %>% 
  ungroup() %>% 
  select(patient,sample,VAR,VAF,sumVAF,group) %>% 
  mutate(status=ifelse(sumVAF>90,"100%","not100%")) %>% 
  arrange(status,group) %>% 
  mutate(number=row_number())


# HC10 (1-3 mutations only)
hc10_vaf_comp <- df %>% filter(patient=="HC10") %>% 
  filter(type=="clonal_normal") %>% 
  filter(vartype!="fe") %>% 
  filter(!(VAR %in% c(freq_recur_var_proj %>% filter(project=="Line1") %>% pull(VAR)))) %>% 
  mutate(across(VAF, ~replace(.,is.na(.),0))) %>% 
  select(patient,VAR,sample,VAF)%>% 
  group_by(sample) %>% 
  mutate(maxVAF=max(VAF)) %>% 
  filter(maxVAF+VAF < 105 | maxVAF==VAF) %>% 
  filter(!(maxVAF>95 & maxVAF!=VAF)) %>% 
  arrange(sample,desc(VAF)) %>% 
  dplyr::slice(1:5)  %>% 
  mutate(secondVAF=sort(VAF,decreasing=TRUE)[2]) %>% 
  mutate(thirdVAF=sort(VAF,decreasing=TRUE)[3]) %>% 
  mutate(maxsecond=maxVAF+secondVAF) %>% 
  mutate(maxthird=maxVAF+secondVAF+thirdVAF) %>% 
  mutate(group=ifelse(maxsecond>105 | is.na(maxsecond), "1mut", ifelse(maxthird>105 | is.na(maxthird), "2mut","3mut"))) %>%
  mutate(sumVAF=ifelse(group=="1mut",maxVAF, ifelse(group=="2mut",maxsecond, maxthird))) %>% 
  dplyr::slice(1:as.numeric(gsub("mut","",group))) %>% 
  ungroup() %>% 
  select(patient,sample,VAR,VAF,sumVAF,group) %>% 
  mutate(status=ifelse(sumVAF>90,"100%","not100%")) %>% 
  arrange(status,group) %>% 
  mutate(number=row_number())


## graph

# tree
tree <- read.tree("~/project/11_Clone_MT/12_Heatmap/02_DB/DB8_Lineage_count_table.txt.nwk")
ggtree(tree)  + geom_rootedge(1) + theme_tree2() + coord_cartesian(xlim=c(0,20))
ggtree(tree)  + geom_rootedge(1) + theme_tree2() + coord_cartesian(xlim=c(1000,10000))

meta_dt <- read_tsv("/home/users/anjisong/project/11_Clone_MT/12_Heatmap/02_DB/Summary_per_sample_210330.txt") %>%
  filter(deadbody=="DB8" & current_final_for_lineage == 'Y') %>% 
  select(lineage_id, sample_id) 

grid.newpage()
pushViewport(viewport(x=0, y=0, width=1, height=1, just=c('left','bottom')))

#### VAF barplot

# to determine the column order of heatmap (sample)
sample_order <- fortify(tree) %>% merge(meta_dt %>% dplyr::rename("label"=lineage_id), by=c("label"), all.x = TRUE) %>% as_tibble() %>% filter(!is.na(sample_id)) %>% arrange(y) %>% pull(sample_id)


pushViewport(viewport(x=0,y=0.35, height=0.15, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))
grid.rect(gp=gpar(lty="solid", fill=FALSE))

grid.yaxis(at=c(0,0.5,1), label=c(0,0.5,1), gp=gpar(fontsize=7))
grid.text('VAF', x=-0.05, rot=90, gp=gpar(fontsize=9))

gt <- db8_vaf_comp %>% 
  group_by(sample) %>% 
  arrange(sample,desc(VAF)) %>% 
  dplyr::slice(1:3) %>% 
  mutate(row=row_number()) %>% 
  arrange(factor(sample, levels=sample_order)) %>% 
  ggplot() + 
  geom_bar(aes(x=factor(sample,levels=sample_order), y=VAF, fill=factor(as.character(row),levels=c("3","2","1"))),stat="identity",position=position_stack()) +
  theme_void() + 
  theme(legend.position="none") + 
  scale_fill_manual(values=c("red","yellow","gray")) + 
  coord_cartesian(ylim=c(5,95))
g <- ggplotGrob(gt)
grid.draw(g)

popViewport(2)



### Fig4d. Svaf vs. turnover --------------------------------------------

normal_mut %>% 
  group_by(project,patient,age,type,tissue2) %>%
  summarise(mean_sumvaf=mean(sumVAF), sd_sumvaf=sd(sumVAF),n=n()) %>% 
  mutate(q1=mean_sumvaf-1.96*sd_sumvaf/sqrt(n), q3=mean_sumvaf+1.96*sd_sumvaf/sqrt(n)) %>% 
  mutate(genyr=ifelse(tissue2=="fibroblast",21, ifelse(tissue2=="colon",14,18))) %>% 
  mutate(total_turnover = genyr*age) %>% 
  ggplot(aes(x=total_turnover,y=mean_sumvaf))  + 
  geom_errorbar(aes(ymin=q1, ymax=q3),col="gray", width=25) +
  geom_point(aes(col=tissue2),size=1.5) + 
  geom_smooth(method='lm',level=0.95, fullrange=TRUE, color="gray30", fill="lightgray") + 
  theme_classic() + 
  scale_color_manual(values=pal_mt) + 
  ggtitle("average VAF sum per patient") 

temp_cor <- normal_mut %>% 
  filter(age!=0) %>% 
  group_by(project,patient,age,type,tissue2) %>%
  summarise(mean_sumvaf=mean(sumVAF), sd_sumvaf=sd(sumVAF),n=n()) %>% 
  mutate(q1=mean_sumvaf-1.96*sd_sumvaf/sqrt(n), q3=mean_sumvaf+1.96*sd_sumvaf/sqrt(n)) %>% 
  mutate(genyr=ifelse(tissue2=="fibroblast",21, ifelse(tissue2=="colon",14,18))) %>% 
  mutate(total_turnover = genyr*age)


# cor.test (turnover-Svaf)
cor.test(temp_cor %>% pull(total_turnover), temp_cor %>% pull(mean_sumvaf))
# cor = 0.787, pva= 2.5e-7



### Fig4e. UV expansion ---------------------------------------------------

uv_df <- df %>% 
  filter(project=="DB") %>% 
  filter(type=="clonal_normal") %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(vartype!="fe") %>%  
  filter(vartype!="fe_gray") %>% 
  mutate(across(VAF, ~replace(., is.na(.), 0))) %>% 
  group_by(patient,sample,age) %>% 
  summarise(max=max(VAF, na.rm=TRUE)) %>% 
  merge(mutsig_df, by=c("patient","sample"), all.x=TRUE) %>% 
  as_tibble()

uv_df %>% 
  mutate(group=ifelse(SBS7num>=10000, "high UV", ifelse(SBS7num < 1500, "low UV", ""))) %>% 
  filter(group!="") %>% filter(!(patient %in% c("2","5","8","9"))) %>% 
  mutate(group_vaf=ifelse(max>=90, "VAF>90", ifelse(max>=75, "VAF>75", ifelse(max>=50, "VAF>50", ifelse(max>=25, "VAF>25", ifelse(max>0,"VAF<25", "VAF<=0")))))) %>% 
  group_by(patient,group,group_vaf) %>% 
  tally() %>% group_by(patient,group) %>% 
  mutate(total=sum(n)) %>% mutate(ratio=n/total*100) %>% 
  mutate(pos=cumsum(ratio)-(0.5*ratio)) %>% 
  ggplot() + 
  geom_bar(aes(fill=group_vaf, x=group, y=n),position="fill", stat="identity") + 
  geom_text(aes(x=group, y=1-pos/100,label=round(ratio))) + 
  theme_light() + ylab("") + 
  facet_grid(patient~.) + 
  scale_fill_brewer(palette="Oranges") + 
  coord_flip() + 
  theme(strip.text.y = element_text(color="black",face="bold"), axis.title.y  = element_blank())



### Fig4f. mutation rate estimation ---------------------------

  ## mitotic turnover 
  load("~/R/R_script/11_Clone_MT/simul_result_mrate_mod2_parse_240130.RData")

  mrate_dist_merge_240129_mod2_v2_uniq_mod <- mrate_dist_merge_240129_mod2_v2_uniq %>% 
    mutate(across(age,~replace(.,is.na(.),0)), across(tissue,~replace(.,is.na(.),"colon"))) %>% 
    arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% ungroup() %>% 
    mutate(num=row_number()) 
  
  # graph
  pat_pos <- seq(1,nrow(mrate_dist_merge_240129_mod2_v2_uniq_mod))
  pat_name <- mrate_dist_merge_240128_mod1_v2_uniq_mod %>% mutate(group=paste0(patient,"(",age,")")) %>% pull(group)
  xaxis_plot_pos = c(0,pat_pos/nrow(mrate_dist_merge_240129_mod2_v2_uniq_mod))
  
  grid.newpage()
  pushViewport(viewport(x=0.1, y=0.1, width=0.8, height=0.5, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))
  
  # background
  for (i in 1:length(pat_pos)){
    if(i %% 2 == 0){bgcolor="gray93"}
    else{bgcolor = FALSE}
    grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
    grid.text(x=xaxis_plot_pos[i]+(1/2)*(1/nrow(mrate_dist_merge_240129_mod2_v2_uniq_mod)), y=-0.05, label = pat_name[i], gp=gpar(fontsize=7), rot=20)
  }
  
  grid.rect(gp=gpar(lty="solid", fill=FALSE))
  grid.yaxis(at=c(0,0.25,0.5,0.75,1), label=c(0,1e-7,2e-7,3e-7,4e-7), gp=gpar(fontsize=7))
  grid.text('mutation rate/bp/generation', x=-0.08, rot=90, gp=gpar(fontsize=9))
  grid.text('Individuals', y=-0.13, gp=gpar(fontsize=9))
  
  grid.segments(x0=0,x1=1, y0=c(0.25,0.5,0.75,1), y1=c(0.25,0.5,0.75,1), gp=gpar(color="gray85", fill=FALSE, lwd=0.2))
  grid.segments(x0=xaxis_plot_pos[c(21,28)], x1=xaxis_plot_pos[c(21,28)],y0=0,y1=1, gp=gpar(color="gray60", fill=FALSE, lwd=0.5))
  
  gt <- mrate_dist_merge_240129_mod2_v2_uniq_mod %>% 
    ggplot(aes(x=reorder(patient,num),y=mrate_avg)) + 
    geom_errorbar(aes(ymin=mrate_min,ymax=mrate_max),col="gray50",width=0.15) +
    geom_point(aes(col=tissue),size=2)  + 
    coord_cartesian(ylim=c(1e-8,3.8e-7)) + 
    theme_void() + theme(legend.position="none")+
    scale_color_manual(values = pal_mt)
  g <- ggplotGrob(gt)
  grid.draw(g)
  
  popViewport(2)
  

  ## homeostatic turnover
  load("~/R/R_script/11_Clone_MT/simul_result_mrate_mod1_parse_240129.RData")

  mrate_dist_merge_240128_mod1_v2_uniq_mod <- mrate_dist_merge_240128_mod1_v2_uniq %>% 
    mutate(across(age,~replace(.,is.na(.),0)), across(tissue,~replace(.,is.na(.),"colon"))) %>% 
    arrange(factor(tissue,levels=c("colon","fibroblast","blood")),age) %>% ungroup() %>% 
    mutate(num=row_number()) 

  
  # graph
  pat_pos <- seq(1,nrow(mrate_dist_merge_240128_mod1_v2_uniq_mod))
  pat_name <- mrate_dist_merge_240128_mod1_v2_uniq_mod %>% mutate(group=paste0(patient,"(",age,")")) %>% pull(group)
  xaxis_plot_pos = c(0,pat_pos/nrow(mrate_dist_merge_240128_mod1_v2_uniq_mod))
  
  grid.newpage()
  pushViewport(viewport(x=0.1, y=0.1, width=0.8, height=0.5, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))
  
  # background
  for (i in 1:length(pat_pos)){
    if(i %% 2 == 0){bgcolor="gray93"}
    else{bgcolor = FALSE}
    grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
    grid.text(x=xaxis_plot_pos[i]+(1/2)*(1/nrow(mrate_dist_merge_240128_mod1_v2_uniq_mod)), y=-0.05, label = pat_name[i], gp=gpar(fontsize=7), rot=20)
  }
  
  grid.rect(gp=gpar(lty="solid", fill=FALSE))
  grid.yaxis(at=c(0,0.2,0.4,0.6,0.8,1), label=c(0,1e-7,2e-7,3e-7,4e-7,5e-7), gp=gpar(fontsize=7))
  grid.text('mutation rate/bp/generation', x=-0.08, rot=90, gp=gpar(fontsize=9))
  grid.text('Individuals', y=-0.13, gp=gpar(fontsize=9))
  
  grid.segments(x0=0,x1=1, y0=c(0.2,0.4,0.6,0.8,1), y1=c(0.2,0.4,0.6,0.8,1), gp=gpar(color="gray85", fill=FALSE, lwd=0.2))
  grid.segments(x0=xaxis_plot_pos[c(21,28)], x1=xaxis_plot_pos[c(21,28)],y0=0,y1=1, gp=gpar(color="gray60", fill=FALSE, lwd=0.5))
  
  gt <- mrate_dist_merge_240128_mod1_v2_uniq_mod %>% 
    ggplot(aes(x=reorder(patient,num),y=mrate_avg)) + 
    geom_errorbar(aes(ymin=mrate_min,ymax=mrate_max),col="gray50",width=0.15) +
    geom_point(aes(col=tissue),size=2)  + 
    coord_cartesian(ylim=c(2e-8,4.8e-7)) + 
    theme_void() + theme(legend.position="none")+
    scale_color_manual(values = pal_mt)
  g <- ggplotGrob(gt)
  grid.draw(g)
  
  popViewport(2)
