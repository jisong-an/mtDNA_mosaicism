################################
#
#   Figure 6 -mtDNA copy number
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
library(gggenes)


### Fig6a. mtCN in normal clone ------------------------------------

sample_info_final_v3 %>% filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(type=="clonal_normal") %>% 
  group_by(tissue2) %>% 
  mutate(med=median(mtCN), mean=mean(mtCN), n=n())%>% 
  ungroup() %>% 
  ggplot() + 
  geom_point(aes(x=reorder(sample,mtCN),y=mtCN), size=1, col="black") + 
  geom_segment(aes(x=10,y=med,xend=n-10,yend=med), color="tomato",size=0.5) +
  facet_grid(. ~ tissue2, scales="free_x") + 
  theme_classic() + 
  theme(legend.position = "none", axis.text.x = element_blank(), strip.background = element_blank(), strip.text = element_text(size=10, colour="black"), axis.ticks.x = element_blank(), axis.title.x = element_blank(), panel.spacing = unit(0.2,"lines")) 

mtcn_df <- sample_info_final_v3 %>% 
  filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(type=="clonal_normal") %>% 
  group_by(patient) %>%
  mutate(mtcn_avg=mean(mtCN)) %>% 
  arrange(factor(tissue2, levels=c("colon","fibroblast","blood")),age,patient,mtCN) %>% 
  ungroup()


# graph
tot_pat <- nrow(mtcn_df)
pat_pos <- mtcn_df %>% 
  group_by(patient) %>% mutate(n=n()) %>% 
  ungroup() %>% select(patient,n) %>% unique() %>% pull(n)
pat_age <- mtcn_df %>% select(patient,age) %>% unique() %>% pull(age)
pat_name <- mtcn_df %>% select(patient) %>% unique() %>% pull(patient)

#xaxis_pos <- c(0)
#for (i in 1:length(pat_pos)){xaxis_pos[i+1]=sum(pat_pos[1:i])}
#xaxis_plot_pos = xaxis_pos/sum(pat_pos)
xaxis_plot_pos = seq(0,length(pat_pos))/length(pat_pos)

grid.newpage()

pushViewport(viewport(x=0,y=0, width=1,height=1,just=c('left','bottom')))

pushViewport(viewport(x=0, y=0.15, width=1, height=0.5, just=c('left','bottom')))
pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))
grid.rect(gp=gpar(lty="solid"))

# background
for (i in 1:length(pat_pos)){
  if(i %% 2 == 0){bgcolor="gray92"}
  else{bgcolor = FALSE}
  grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
}

grid.yaxis(at=c(0,0.2,0.4,0.6,0.8,1), label=c(0,500,1000,1500,2000,2500), gp=gpar(fontsize=7))
grid.xaxis(label = FALSE, at=c(0,1))

grid.text('mtDNA copy number', x=-0.03, rot=90, gp=gpar(fontsize=9))
grid.text('clones', y=-0.1, gp=gpar(fontsize=9))


# y-value : mtCN
yrange=2500
yaxis_plot_pos = seq(0,1,by=1/yrange)
y_value <- mtcn_df %>% mutate(n=mtCN/yrange) %>% 
  arrange(factor(tissue2, levels=c("colon","fibroblast","blood")),age,patient,mtCN) %>% 
  pull(n)  #mtcn
y_value2 <- mtcn_df %>% select(patient,mtcn_avg) %>% unique() %>% mutate(avg=mtcn_avg/yrange) %>% 
  arrange(factor(tissue2, levels=c("colon","fibroblast","blood")),age,patient,mtCN) %>% 
  pull(avg)  #mtcn_avg

# x : each clone
xaxis_data <- mtcn_df %>% mutate(order=row_number()/tot_pat) %>% pull(order)


# dot plot
now=0
for (i in 1:length(pat_pos)){
  clone_n=pat_pos[i]
  xaxis_sub_pos = seq(xaxis_plot_pos[i],xaxis_plot_pos[i+1],length.out=clone_n+4)
  for (j in 1:clone_n){
    grid.points(x=xaxis_sub_pos[j+2],y=y_value[now+j], pch=16, size=unit(0.4,'char'), gp=gpar(col="black"))  # mtCN
  }
  now = now+clone_n
  grid.segments(x0=xaxis_sub_pos[5],x1=xaxis_sub_pos[clone_n-5], y0=y_value2[i], y1=y_value2[i], gp=gpar(col="red", fill=FALSE, lwd=1))
  grid.text(x = (xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=-0.04, label = pat_name[i], gp=gpar(fontsize=7), rot=30)
  
}

popViewport(3)



### Fig 6b. large SV in normal clone ---------------------------------

load("~/R/R_script/11_Clone_MT/RNAseq_HC21_result.RData")

hc21_count_mod <- hc21_count %>% 
  mutate(gene_group=ifelse(grepl("^MT-",gene_name), gene_name, "n")) %>% 
  group_by(sample,gene_group, length) %>% 
  summarise(sum=sum(expected_count)) %>% 
  mutate(group=ifelse(sample%in%group1, "group1", "group2")) %>% 
  filter(group=="group2") %>% 
  ungroup() %>% 
  group_by(sample) %>% 
  mutate(ratio=sum/sum(sum)) %>% 
  filter(gene_group != "n") %>% 
  mutate(group=ifelse(sample=="HC21-16","largeSV", group)) %>% 
  mutate(ratio_mod=ratio/length)


depth_hc21 <- read_tsv("~/project/11_Clone_MT/02_Pileup/01_Line1/HC21-16.MT.all.tsv")

grid.newpage()
pushViewport(viewport(x=0,y=0, width=1, height=1, just=c('left','bottom')))
pushViewport(viewport(x=0.1,y=0.5, width=0.8, height=0.4, just=c('left','bottom')))
pushViewport(viewport(x=0.5,y=0.5, width=0.9, height=0.9))

grid.yaxis(at=c(0,0.33,0.66,1), label=c(0,5000,10000,15000), gp=gpar(fontsize=7))
grid.text('mtDNA read depth', x=-0.05, rot=90, gp=gpar(fontsize=9))

for (i in 1:nrow(depth_hc21)){
  grid.points(x=depth_hc21[[i,"POS"]]/16569, y=depth_hc21[[i,"Depth"]]/15000, pch=16, size=unit(0.25, 'char'))
}

popViewport(3)


depth_hc21 %>% ggplot() + geom_line(aes(x=POS,y=Depth, group=1)) + coord_cartesian(ylim=c(0,15000)) + theme_classic()


# HC06 deletion
depth_hc06 <- read_tsv("~/project/11_Clone_MT/02_Pileup/01_Line1/HC06-14.MT.all.tsv")
depth_hc06 %>% ggplot() + geom_line(aes(x=POS,y=Depth, group=1)) + coord_cartesian(ylim=c(0,10000)) + theme_classic()


# genes
mito_genes <- as_tibble(read.table("/home/users/kimin/projects/12_MT/00_reference/mito_genes.bed", col.names = c("chrom","start","end","gene","dummy","strand","type"), colClasses = c("character","integer","integer","character","integer","character","character")))
mito_genes$direction <- ifelse(mito_genes$strand == "+", 1, -1)
mito_genes$gene_label <- ifelse((mito_genes$end-mito_genes$start) < 350, "", mito_genes$gene)

gene_color_vector <- pal_jama()(7)
names(gene_color_vector) <- sort(unique(mito_genes$type))

genes <- ggplot(mito_genes, aes(xmin=start, xmax=end, y=chrom, fill=type, label=gene, forward=direction)) + 
  geom_gene_arrow(arrow_body_height = unit(5,"mm"), arrowhead_height=unit(5, "mm"), arrowhead_width=unit(1.5, "mm")) +
  geom_gene_label(align = 'centre') +
  facet_wrap(~ chrom, scales='free', ncol=1) +
  scale_fill_jama() +
  theme_genes() + guides(fill="none") + 
  theme(legend.position='bottom') + ylab('') + xlab('') +
  scale_x_continuous(limits = c(0,16569), expand = c(0, 0), breaks = seq(0, 16000, by = 4000))



### Fig6c. colon cancer & normal mutation (count) ------------------------------

colon_type_mut <- df %>% 
  filter(tissue2=="colon") %>% 
  filter(vartype!="fe" & vartype!="fe_gray") %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(type!="bulk_normal") %>% 
  filter(!is.na(VAR)) %>% 
  filter(project %in% c("Abortus","Line1")) %>% 
  group_by(project,patient,age,sample,type) %>% 
  mutate(sumVAF=sum(VAF), n=n(), maxVAF=max(VAF)) %>% 
  group_by(project,patient,age,sample,type,sumVAF,maxVAF) %>% 
  tally() %>% 
  merge(sample_info_final_v3 %>% filter(type!="bulk_normal")  %>% 
          filter(project %in% c("Abortus","Line1")) %>% 
          dplyr::select(project,patient,sample,age,type),by=c("project","patient","sample","age","type"),all=TRUE) %>%
  as_tibble() %>% 
  mutate(across(n,~replace(.,is.na(.),0)), across(sumVAF,~replace(.,is.na(.),0)), across(maxVAF, ~replace(.,is.na(.),0)))


colon_type_mut %>%  
  group_by(project,patient,age,type) %>%
  summarise(mean_n=mean(n), sd_n=sd(n),n=n()) %>% 
  mutate(q1=mean_n-1.96*sd_n/sqrt(n), q3=mean_n+1.96*sd_n/sqrt(n)) %>% 
  mutate(genyr=14.3) %>% 
  mutate(total_turnover = genyr*age) %>% 
  ggplot(aes(x=total_turnover,y=mean_n,col=type)) + 
  geom_errorbar(aes(ymin=q1,ymax=q3),col="gray",width=30)+
  geom_point() + 
  geom_smooth(method='lm', se=FALSE) + 
  theme_classic() + 
  ggtitle("average mutation count per patient") + 
  scale_color_manual(values=c(pal_mt[[1]],"#B83535")) + 
  scale_x_continuous(breaks=seq(0,1200,400))


## paired statistics
temp <- colon_type_mut %>% 
  filter(project!="CRC") %>% filter(project!="MUTYH") %>% 
  group_by(project,patient,age,type) %>%
  summarise(mean_n=mean(n), sd_n=sd(n),n=n()) %>% 
  mutate(q1=mean_n-1.96*sd_n/sqrt(n), q3=mean_n+1.96*sd_n/sqrt(n)) %>% 
  mutate(genyr=14.3) %>% 
  mutate(total_turnover = genyr*age) %>% filter(age!=0)  


wilcox.test(temp %>% filter(type=="tumor") %>% pull(mean_n), 
            temp %>% filter(type=="clonal_normal") %>% pull(mean_n), paired=TRUE, alternative="greater")
# p=0.0301



### Fig6d. colon cancer & normal mutation (Svaf) ------------------------------

colon_type_mut <- df %>% 
  filter(tissue2=="colon") %>% 
  filter(vartype!="fe" & vartype!="fe_gray") %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(type!="bulk_normal") %>% 
  filter(!is.na(VAR)) %>% 
  filter(project %in% c("Abortus","Line1")) %>% 
  group_by(project,patient,age,sample,type) %>% 
  mutate(sumVAF=sum(VAF), n=n(), maxVAF=max(VAF)) %>% 
  group_by(project,patient,age,sample,type,sumVAF,maxVAF) %>% 
  tally() %>% 
  merge(sample_info_final_v3 %>% filter(type!="bulk_normal")  %>% 
          filter(project %in% c("Abortus","Line1")) %>% 
          dplyr::select(project,patient,sample,age,type),by=c("project","patient","sample","age","type"),all=TRUE) %>%
  as_tibble() %>% 
  mutate(across(n,~replace(.,is.na(.),0)), across(sumVAF,~replace(.,is.na(.),0)), across(maxVAF, ~replace(.,is.na(.),0)))


# mean Svaf
colon_type_mut %>% 
  group_by(project,patient,age,type) %>%
  summarise(mean_sumvaf=mean(sumVAF), sd_sumvaf=sd(sumVAF),n=n()) %>% 
  mutate(q1=mean_sumvaf-1.96*sd_sumvaf/sqrt(n), q3=mean_sumvaf+1.96*sd_sumvaf/sqrt(n)) %>% 
  mutate(genyr=14.3) %>% 
  mutate(total_turnover = genyr*age) %>% 
  ggplot(aes(x=total_turnover,y=mean_sumvaf,col=type)) + 
  geom_errorbar(aes(ymin=q1,ymax=q3),col="gray",width=30)+
  geom_point() + 
  geom_smooth(method='lm', se=FALSE) + 
  theme_classic() + 
  ggtitle("average VAF sum per patient") + 
  scale_color_manual(values=c(pal_mt[[1]],"#B83535"))+
  scale_x_continuous(breaks=seq(0,1200,400))


## paired statistics
temp <- colon_type_mut %>% 
  filter(project!="CRC") %>% filter(project!="MUTYH") %>% 
  group_by(project,patient,age,type) %>%
  summarise(mean_sumvaf=mean(sumVAF), sd_sumvaf=sd(sumVAF),n=n()) %>% 
  mutate(q1=mean_sumvaf-1.96*sd_sumvaf/sqrt(n), q3=mean_sumvaf+1.96*sd_sumvaf/sqrt(n)) %>% 
  mutate(genyr=14.3) %>% 
  mutate(total_turnover = genyr*age)  %>% filter(age!=0)  


wilcox.test(temp %>% filter(type=="tumor") %>% pull(mean_sumvaf), 
            temp %>% filter(type=="clonal_normal") %>% pull(mean_sumvaf), paired=TRUE, alternative="greater")
# p=0.0008469



### Fig6e. truncating mutation between turmo & normal ----------------------

  df %>% 
    filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
    filter(consequence=="truncating") %>% 
    filter(VAF>60) %>% 
    select(VAR,VAF,type) %>% 
    rbind(pcawg_var_df %>% filter(consequence=="truncating") %>% 
            filter(VAF>60) %>% select(VAR,VAF) %>% mutate(type="tumor")) %>% 
    as_tibble() %>% 
    group_by(type) %>% tally() %>% 
    merge(mut_count_type_pcawg %>% dplyr::rename("mutN"=n),by=c("type")) %>% 
    mutate(ratio=n/mutN) %>% 
    ggplot() + 
    geom_bar(aes(x=type,y=ratio,fill=type),stat="identity") + 
    geom_text(aes(x=type,y=ratio,label=paste0(n,"/",mutN)),vjust=-0.5) + 
    theme_classic() + 
    ggtitle("truncating mutation (VAF>60%)")
  
  
  # fisher test
  count <- c(5677-15,344-7,15,7)
  group <- c("N","C","N","C")  # N : normal, C : cancer
  group2 <- c("T","T","X","X") # T : total var, X : truncating
  dat <- data.frame(group,group2,count)
  tab <- xtabs(count~group+group2,data=dat)
  fisher.test(tab)
  # 0.000151
  # 0.0203 vs 0.0026



### Fig6f. colon tumor & normal mtCN ---------------------------------------

mtcn_df_diploid <- sample_info_final_v3 %>% 
  filter(project=="Line1") %>% 
  filter(type!="bulk_normal") %>% 
  mutate(mtCN=mtDNA_depth/nDNA_depth*2) %>% 
  group_by(patient,type) %>% 
  mutate(mtCN_avg=mean(mtCN))

order_pat <- mtcn_df_diploid %>% filter(type=="clonal_normal") %>% select(patient,mtCN_avg) %>% unique() %>% arrange(mtCN_avg) %>% pull(patient)

ggplot(data=subset(mtcn_df_diploid, type=="clonal_normal"), aes(x=factor(patient,levels=order_pat),y=mtCN, col=type)) + 
  geom_jitter(width=0.1,size=1.25) + 
  geom_boxplot(alpha=0, width=0.5, col="gray30") + 
  geom_point(data=subset(mtcn_df_diploid, type=="tumor"), aes(x=factor(patient,levels=order_pat),y=mtCN, col=type), size=3.5, pch=18)+
  theme_classic() + scale_color_manual(values=c(pal_mt[[1]],"#B83535"))  + 
  theme(axis.text.x=element_text(angle=20)) + xlab("Individuals")



### Fig 6g. tumor cell fraction-tumor mtCN ---------------------------------------------

tumor_cninfo <- read_tsv("~/project/11_Clone_MT/info/01_Line1/tumor_CN_info.tsv") 
mtcn_df_diploid_purity <- mtcn_df_diploid %>% 
  filter(type=="tumor") %>% 
  merge(tumor_cninfo %>% dplyr::rename("sample"=sampleID) %>% 
          select(sample,purity,ploidy),by=c("sample")) %>% 
  as_tibble()

# tumor CN - purity (tumor cell fraction)
mtcn_df_diploid_purity %>% 
  ggplot(aes(x=purity,y=mtCN)) +  
  geom_smooth(method="lm", level=0.95, fullrange=TRUE, color="gray", fill="gray90") + 
  geom_point(size=2,col="#B83535") + 
  theme_classic() + xlim(c(0,1))

cor.test(mtcn_df_diploid_purity$purity, mtcn_df_diploid_purity$mtCN)  # pval = 0.000573, cor=0.715, r=0.715

# univariable linear regression analysis
lm_fit <- lm(mtCN~purity, data=mtcn_df_diploid_purity)
summary(lm_fit)  

# =mtcn = 211.7 + 1053.9*purity
# purity=1; mtcn = 1265.6
