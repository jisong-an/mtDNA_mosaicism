################################
#
#   Figure 2 -HetFE variants
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


### Fig2c. DB2 16400 C>T  ---------------------------------------

tree <- read.tree("~/project/11_Clone_MT/12_Heatmap/02_DB/DB2_Lineage_count_table.txt.nwk")
meta_dt <- read_tsv("/home/users/anjisong/project/11_Clone_MT/12_Heatmap/02_DB/Summary_per_sample_210330.txt") %>% 
  filter(deadbody=="DB2" & current_final_for_lineage == 'Y') %>% 
  select(lineage_id, sample_id)

#### tree graph
tree$edge.length <- log(tree$edge.length+1,2)
tree$root.edge <- 1
p <- ggtree(tree)  + 
  geom_tippoint(size=1) + 
  layout_dendrogram() + geom_rootedge()

p$data$branch_label <- round(2**p$data$branch.length-1,0)
p$data <- p$data %>% mutate(branch_label=ifelse(branch_label >= 30, "", branch_label))
p + geom_text(aes(x=branch, y=y,label=branch_label), hjust=-1, size=2.5) 


#### VAF barplot

# to determine the column order of heatmap (sample)
sample_order <- fortify(tree) %>% merge(meta_dt %>% rename("label"=lineage_id), by=c("label"), all.x = TRUE) %>% as_tibble() %>% filter(!is.na(sample_id)) %>% arrange(y) %>% pull(sample_id)


data_n <- df %>% filter(patient=="2" & POS==16400) %>% 
  select(sample,VAF) %>% 
  merge(sample_info_final_v3 %>% filter(patient=="2") %>% 
          select(sample),by=c("sample"),all=TRUE) %>% as_tibble() %>% 
  mutate(across(VAF,~replace(.,is.na(.),0)))%>% arrange(factor(sample,levels=sample_order))

# y-value : VAF
#yrange = data_n %>% pull(n) %>% max()
yaxis_plot_pos = c(0,0.5,1)
y_value <- data_n %>% mutate(VAF=VAF/100) %>% pull(VAF)

# x : each clone
xaxis_data <- data_n %>% mutate(order=row_number()/nrow(data_n)) %>% pull(order)

### draw graph
grid.newpage()

pushViewport(viewport(x=0, y=0.1, width=1, height=0.15, just=c('left','bottom')))
pushViewport(viewport(x=1, y=0.5, height = 0.9, width=0.9, just=c('right')))
grid.rect(gp=gpar(lty="solid"))

# axis
#grid.xaxis(at=xaxis_plot_pos, label=FALSE)
grid.yaxis(at=yaxis_plot_pos, label=yaxis_plot_pos, gp = gpar(fontsize=8))
grid.text('VAF', x=-0.07, rot=90)


# bar graph
for (i in 1:nrow(data_n)){
  grid.rect(x=xaxis_data[i],y=0, just=c('right','bottom'), height=y_value[i], width=xaxis_data[1], gp=gpar(fill=pal_mt[2], col="lightgray", lwd=0.5))
  grid.text(round(y_value[i],2), x=xaxis_data[i]-xaxis_data[1]/2, y=-0.1, gp=gpar(fontsize=7.5))
}

popViewport(2)



### Fig2d. HC19 7496 T>C ----------------------------------------------------

tree <- read.tree("~/project/11_Clone_MT/12_Heatmap/01_Line1/HC19.length.nwk")
# to determine the column order of heatmap (sample)
sample_order <- fortify(tree) %>% filter(!is.na(label)) %>% arrange(y) %>% pull(label)


#### tree graph
tree$edge.length <- log(tree$edge.length+1,2)
tree$root.edge <- 1
p <- ggtree(tree)  + 
  geom_tippoint(size=1) + 
  layout_dendrogram() + geom_rootedge()

p$data$branch_label <- round(2**p$data$branch.length-1,0)
p$data <- p$data %>% mutate(branch_label=ifelse(branch_label >= 30, "", branch_label))
p + geom_text(aes(x=branch, y=y,label=branch_label), hjust=-1, size=2.5)



data_n <- df %>% filter(patient=="HC19" & POS==7496) %>% filter(type!="bulk_normal") %>%  
  select(sample,VAF) %>% 
  merge(sample_info_final_v3 %>% filter(patient=="HC19") %>% filter(type!="bulk_normal") %>% 
          select(sample),by=c("sample"),all=TRUE) %>% as_tibble() %>% 
  mutate(across(VAF,~replace(.,is.na(.),0))) %>% arrange(factor(sample,levels = sample_order))


# y-value : VAF
#yrange = data_n %>% pull(n) %>% max()
yaxis_plot_pos = c(0,0.5,1)
y_value <- data_n %>% mutate(VAF=VAF/100) %>% pull(VAF)

# x : each clone
xaxis_data <- data_n %>% mutate(order=row_number()/nrow(data_n)) %>% pull(order)

### draw graph
grid.newpage()

pushViewport(viewport(x=0, y=0.1, width=1, height=0.15, just=c('left','bottom')))
pushViewport(viewport(x=1, y=0.5, height = 0.9, width=0.9, just=c('right')))
grid.rect(gp=gpar(lty="solid"))

# axis
#grid.xaxis(at=xaxis_plot_pos, label=FALSE)
grid.yaxis(at=yaxis_plot_pos, label=yaxis_plot_pos, gp = gpar(fontsize=8))
grid.text('VAF', x=-0.07, rot=90)


# bar graph
for (i in 1:nrow(data_n)){
  grid.rect(x=xaxis_data[i],y=0, just=c('right','bottom'), height=y_value[i], width=xaxis_data[1], gp=gpar(fill=pal_mt[1], col="lightgray", lwd=0.5))
  grid.text(round(y_value[i],2), x=xaxis_data[i]-xaxis_data[1]/2, y=-0.1, gp=gpar(fontsize=7.5))
}

popViewport(2)



### Fig2e. blood VAF - pseudo VAF correlation --------------------------------

df_blood_vaf_fe <- df_blood_vaf %>% 
  merge(rbind(fe_blood,uniq_blood), by=c("project","patient","VAR")) %>% 
  as_tibble()  %>% select(-pseudoVAF) %>% 
  merge(df %>% select(patient,VAR,pseudoVAF) %>% unique(),by=c("patient","VAR"),all.x=TRUE) %>% 
  as_tibble()

df_blood_vaf_fe %>% 
  ggplot(aes(x=pseudoVAF/100, y=bloodVAF/100))  + 
  geom_smooth(method="lm", level=0.95, fullrange=TRUE, color="blue", fill="lightblue") + 
  geom_point(size=1.5) + 
  theme_classic() + xlab("pseudo VAF") + ylab("blood VAF")

cor.test(df_blood_vaf_fe %>% pull(bloodVAF), df_blood_vaf_fe %>% pull(pseudoVAF))
# r= 0.967, p = 2.299e-16



### Fig2f. pseudoVAF & clone count -----------------------------------

## pseudoVAF distribution per clone count
pseudoVAF_clone_all <- pseudoVAF_clone %>% 
  merge(df %>% filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
          select(patient,patientN) %>% unique(), by=c("patient","patientN"), all=TRUE) %>% 
  as_tibble()

pseudoVAF_clone_all_n <- pseudoVAF_clone_all %>% 
  mutate(vaf_group=ifelse(log2_pseudoVAF>5, ">32%", 
                          ifelse(log2_pseudoVAF>4, ">16%", 
                                 ifelse(log2_pseudoVAF>3,">8%", 
                                        ifelse(log2_pseudoVAF>2,">4%", 
                                               ifelse(log2_pseudoVAF>1, ">2%", 
                                                      ifelse(log2_pseudoVAF>0, ">1%", 
                                                             ifelse(log2_pseudoVAF>-1,">0.5%",
                                                                    ifelse(log2_pseudoVAF>-2, ">0.25%", 
                                                                           ifelse(log2_pseudoVAF>-3,">0.125%", 
                                                                                  ifelse(log2_pseudoVAF>-4,">0.0625%", 
                                                                                         ifelse(log2_pseudoVAF>-5, ">0.0312%", 
                                                                                                ifelse(log2_pseudoVAF>-6, ">0.0156%",
                                                                                                       ifelse(log2_pseudoVAF>-7, ">0.0088%", ">0.0044%")))))))))))))) %>% 
  group_by(patient, patientN, vaf_group) %>% 
  tally() %>% mutate(n=ifelse(is.na(vaf_group), NA, n))


ggplot(data=pseudoVAF_clone_all_n %>%
         rbind(tibble(patient="", patientN=NA, vaf_group=">16%",n=NA)) %>% 
         mutate(group=paste0(patientN,"\n(",patient,")")), 
       aes(x=reorder(group,patientN), 
           y=factor(vaf_group, 
                    levels=rev(c(">32%",">16%",">8%",">4%",">2%",">1%",">0.5%",">0.25%",">0.125%",">0.0625%",">0.0312%",">0.0156%",">0.0088%",">0.0044%"))))) + 
  geom_point(shape=21, fill="goldenrod1", color="gray90", size=6) +
  geom_text(aes(label=n), size=3) + 
  geom_hline(yintercept=seq(1.5,13.5,1), color="gray90") + 
  geom_vline(xintercept=seq(1.5,30.5,1), color="gray90") + 
  coord_cartesian(xlim=c(1,31),ylim=c(0.8,14)) + 
  theme_classic() + 
  xlab("number of clone") + ylab("pseudo VAF") 

# patient ratio
pseudoVAF_patient <- pseudoVAF_clone_all_n %>% 
  mutate(total_pat=ifelse(vaf_group %in% c(">0.0625%"), 27, 
                          ifelse(vaf_group %in% c(">0.0312%"), 25, 
                                 ifelse(vaf_group %in% c(">0.0156%",">0.0088%"), 5, 
                                        ifelse(vaf_group %in% c(">0.0044%"), 3, 31)))))   %>% 
  group_by(vaf_group, total_pat) %>% 
  mutate(totvar=sum(n)) %>% 
  filter(!is.na(vaf_group)) %>% 
  select(vaf_group, total_pat, totvar) %>% 
  unique() %>% 
  mutate(ratio=totvar/total_pat) %>% 
  arrange(factor(vaf_group, levels=c(">32%",">16%",">8%",">4%",">2%",">1%",">0.5%",">0.25%",">0.125%",">0.0625%",">0.0312%",">0.0156%",">0.0088%",">0.0044%"))) %>% 
  ungroup()

# weighted distribution
# df_autism_n2 from autism_analysis.R
ggplot(data=pseudoVAF_patient, aes(x=factor(vaf_group, 
                                            levels=rev(c(">32%",">16%",">8%",">4%",">2%",">1%",">0.5%",">0.25%",">0.125%",">0.0625%",">0.0312%",">0.0156%",">0.0088%",">0.0044%"))),
                                   y=ratio))  + 
  geom_hline(yintercept=seq(0,2,0.25), color="gray90") + 
  geom_bar(stat="identity", fill="chocolate4") + 
  geom_point(data=df_autism_n2, aes(x=factor(vaf_group, levels=rev(c(">32%",">16%",">8%",">4%",">2%",">1%",">0.5%",">0.25%",">0.125%",">0.0625%",">0.0312%",">0.0156%",">0.0088%",">0.0044%"))), y=ratio), col="gray20", size=2) + 
  coord_flip() + 
  theme_classic() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) + 
  ylim(0,2) + ylab("Avg. mutation count per patient")


save(pseudoVAF_clone, pseudoVAF_clone_all, pseudoVAF_clone_all_n, pseudoVAF_patient, file="~/R/R_script/11_Clone_MT/pseudoVAf_clone_heatmap.RData")


## tissue & clone count bar 
pat_df <- sample_info_final_v3 %>% filter(project %in% c("Abortus","Line1","DB","Hblood","blood_new")) %>% filter(type=="clonal_normal") %>% group_by(patient,tissue2) %>% tally() %>% arrange(n)
pat_tissue <- pat_df %>% pull(tissue2)
pat_clone <- pat_df %>% pull(n)

xaxis_plot_pos = seq(0,31)/31

grid.newpage()
pushViewport(viewport(x=0,y=0, width=1,height=1,just=c('left','bottom')))

# clone count
pushViewport(viewport(x=0,y=0.25, height=0.05, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))

grid.text('# of clone', x=-0.03, gp=gpar(fontsize=8))
col_fun = colorRamp2(c(10,50,400),c("white",pal_jco()(10)[2],pal_jama()(10)[2]))
#col_fun = colorRamp2(c(0,30,100),c("#FDFBF5","#FF97B5","#6E2A85"))

for (i in 1:length(pat_clone)){
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=col_fun(pat_clone[i]), col=FALSE))
  grid.text(pat_clone[i], x=(xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=0.5, gp=gpar(fontsize=6.5))
}
popViewport(3)

## tissue
pushViewport(viewport(x=0,y=0.2, height=0.05, just=c('left','bottom')))
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
popViewport(2)



### Fig2g. VAF in offspring (family bulk tissue) -----------------------------------

df_family_merge_check_filt_all_mod_v2 %>% 
  filter(group!="mother_only") %>% 
  filter(patient!="ST304") %>% 
  rbind(df_family_merge_check_filt_all_mztwin %>% 
          filter(patient=="ST304") %>% ungroup() %>% dplyr::slice(1:2)) %>% 
  mutate(vaf_group=ifelse(VAF>32,">32",ifelse(VAF>16,"16-32",ifelse(VAF>8,"8-16",ifelse(VAF>4,"4-8",ifelse(VAF>2,"2-4",ifelse(VAF>1,"1-2", "0.5-1"))))))) %>% 
  group_by(vaf_group,group) %>% 
  tally() %>% 
  mutate(ratio=n/407) %>% 
  ggplot(aes(x=factor(vaf_group,levels=rev(c(">32","16-32","8-16","4-8","2-4","1-2","0.5-1"))),y=ratio)) + 
  geom_hline(yintercept=seq(0,0.3,0.1),col="lightgray")+
  geom_bar(aes(fill=group),stat="identity",position=position_stack()) +
  xlab("VAF in bulk tissue") + 
  coord_flip() + 
  scale_fill_manual(values=pal_simpsons()(10)[c(2,1)]) + 
  theme_classic()
