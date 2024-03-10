################################
#
#   Figure 5 -functional impact
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


### Fig5a. dN/dS ------------------------------------

# obs_merge_consq from result_simul_selection_230724.R
obs_merge_consq_mod <- obs_merge_consq %>% 
  ungroup() %>% 
  filter((nonsynonymous_n > 10 & consequence=="nonsynonymous") | (nonsense_n > 10 & consequence=="nonsense")) %>% 
  filter(context!="A>C") %>% 
  select(tissue,consequence,context,dnds,everything()) %>% 
  group_by(tissue,consequence) %>% 
  summarise(dnds_avg=mean(dnds))  # merge all context

simul_result_merge_consq_mod <- simul_result_merge_consq %>% 
  group_by(consequence,tissue) %>% 
  summarise(q1_avg=mean(q1),q2_avg=mean(q2))


# graph
obs_merge_consq_mod %>% 
  merge(simul_result_merge_consq_mod, by=c("consequence","tissue")) %>% 
  as_tibble() %>% 
  mutate(tissue=ifelse(tissue=="fibr","fibroblast",tissue)) %>% 
  ggplot() + 
  geom_errorbar(aes(x=factor(tissue,levels=c("colon","fibroblast","blood")), ymin=q1_avg, ymax=q2_avg), width=0.1, color="gray40")+
  geom_point(aes(x=factor(tissue,levels=c("colon","fibroblast","blood")),y=dnds_avg, col=tissue),size=4,shape=18)+
  theme_classic() + 
  facet_wrap(~factor(consequence,levels=c("nonsynonymous","nonsense"))) +
  scale_color_manual(values=pal_mt)+
  ylab("dN/dS") + xlab("tissue") + coord_cartesian(ylim=c(0,2))+
  theme(strip.background = element_rect(fill="lightgray", color=NA), 
        strip.text.x = element_text(face="bold"), title = element_text(size=10))




### Fig5b. VAF & consequence -----------------------------------------

df %>% 
  filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(!(vartype %in% c("fe","fe_gray"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(type=="clonal_normal") %>% 
  filter(!is.na(VAR))  %>% 
  mutate(consequence = ifelse(Consequence=="downstream_gene_variant" | Consequence=="upstream_gene_variant","intergenic",ifelse(Consequence=="missense_variant","missense",ifelse(Consequence=="synonymous_variant","synonymous",ifelse(Consequence=="inframe_deletion" | Consequence=="inframe_insertion","inframe indel",ifelse(Consequence=="stop_gained" | Consequence=="frameshift_variant","truncating", ifelse(Consequence=="stop_lost","nonstop",Consequence))))))) %>% 
  mutate(group_vaf=ifelse(VAF>90, "VAF>90", ifelse(VAF>75, "VAF>75", ifelse(VAF>50, "VAF>50",ifelse(VAF>25, "VAF>25", ifelse(VAF>10, "VAF>10", "VAF<10")))))) %>% 
  filter(consequence %in% c("intergenic","missense","synonymous","truncating")) %>% 
  filter(consequence!="intergenic") %>% 
  group_by(consequence, group_vaf, tissue2) %>% 
  tally() %>% 
  group_by(consequence, tissue2) %>% 
  mutate(tot=sum(n)) %>% 
  mutate(ratio=round(n/tot,2))%>% 
  ungroup() %>% 
  filter(group_vaf!="VAF<10") %>% 
  group_by(consequence,tissue2) %>% 
  ggplot() + 
  geom_bar(aes(fill=factor(tissue2,levels=c("colon","fibroblast","blood")), x=factor(consequence,levels=c("synonymous","missense","truncating")), y=ratio, group=group_vaf, alpha=group_vaf), position=position_stack(), stat="identity") + 
  geom_text(aes(x=factor(consequence,levels=c("synonymous","missense","truncating")), y=ratio,label=round(ratio,2), group=group_vaf),size=3, position=position_stack(0.5)) + 
  facet_wrap(~factor(tissue2, levels=c("colon","fibroblast","blood"))) + 
  theme_light() + ylab("") +  
  scale_fill_manual(values=c(brewer.pal("Blues",n=6)[6],brewer.pal("Oranges",n=6)[6],brewer.pal("Greens",n=6)[6])) + 
  scale_alpha_manual(values=c(0.1,0.3,0.5,0.7,1)) + 
  theme(strip.text.x = element_text(color="black"),panel.grid = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle=20))




### Fig5c. Truncating & TPM (HC13-10,13) ----------------------------------

# truncating - HC13-10, 13
rna_count_hc13_trunc <- rna_count_df_mod_hc13_trunc %>% as.data.frame() %>% 
  rownames_to_column(var="gene") %>% as_tibble() %>% 
  gather(2:ncol(.),key="sample",value="count") %>% 
  group_by(sample) %>% mutate(sum=sum(count)) %>% mutate(ratio=count/sum) %>% 
  filter(grepl("^MT-",gene)) %>% 
  mutate(group=ifelse(sample %in% c("HC01-05","HC04-05","HC05-16","HC13-10","HC13-13","HC17-01"),"truncating", ifelse(sample %in% c("HC13-12"),"missense in CO3","normal")))
#filter(!(gene%in%c("MT-ND6","MT-ATP8"))) 

rna_count_hc13_trunc_ratio <- rna_count_hc13_trunc %>% 
  merge(rna_count_hc13_trunc %>% filter(group!="truncating") %>% 
          filter(group!="missense in CO3") %>% group_by(gene) %>% mutate(mean=mean(ratio)) %>% 
          select(gene,mean) %>% unique(),by=c("gene")) %>% 
  as_tibble() %>% 
  mutate(FC=ratio/mean) %>% 
  mutate(log_FC=log2(FC))

rna_count_hc13_trunc_ratio %>% 
  filter(group=="normal") %>% 
  ggplot(aes(x=gene,y=log_FC,col=group)) + 
  geom_hline(yintercept=c(-0.5,0,0.5,1,1.5,2), col="gray90") + 
  geom_jitter(width=0.1) + 
  geom_boxplot(alpha=0, col="gray30") + 
  geom_point(data=subset(rna_count_hc13_trunc_ratio, group=="truncating"), aes(x=gene,y=log_FC), size=3.5, pch=18) + 
  geom_point(data=subset(rna_count_hc13_trunc_ratio, group=="missense in CO3"), aes(x=gene,y=log_FC), size=3.5, pch=18) + 
  scale_color_cosmic(palette="signature_substitutions") + 
  theme_classic() + coord_cartesian(ylim=c(-0.75,2))




### Fig5d. Truncating & TPM (HC17-01) -------------------------------

# truncating - HC17-01
rna_count_hc17_trunc <- rna_count_df_mod_hc17_trunc %>% as.data.frame() %>% 
  rownames_to_column(var="gene") %>% as_tibble() %>% 
  gather(2:ncol(.),key="sample",value="count") %>% 
  group_by(sample) %>% mutate(sum=sum(count)) %>% mutate(ratio=count/sum) %>% 
  filter(grepl("^MT-",gene)) %>% 
  mutate(group=ifelse(sample %in% c("HC01-05","HC04-05","HC05-16","HC13-10","HC13-13","HC17-01"),"truncating", ifelse(sample %in% c("HC17-07"),"missense in ND4","normal")))
#filter(!(gene%in%c("MT-ND6","MT-ATP8"))) 

rna_count_hc17_trunc_ratio <- rna_count_hc17_trunc %>% 
  merge(rna_count_hc17_trunc %>% filter(group!="truncating") %>% 
          filter(group!="missense in ND4") %>% group_by(gene) %>% mutate(mean=mean(ratio)) %>% 
          select(gene,mean) %>% unique(),by=c("gene")) %>% 
  as_tibble() %>% 
  mutate(FC=ratio/mean) %>% 
  mutate(log_FC=log2(FC))

rna_count_hc17_trunc_ratio %>% 
  filter(group=="normal") %>% 
  ggplot(aes(x=gene,y=log_FC,col=group)) + 
  geom_hline(yintercept=c(-0.5,0,0.5,1), col="gray90") + 
  geom_jitter(width=0.1) + 
  geom_boxplot(alpha=0, col="gray30") + 
  geom_point(data=subset(rna_count_hc17_trunc_ratio, group=="truncating"), aes(x=gene,y=log_FC), size=3.5, pch=18) + 
  geom_point(data=subset(rna_count_hc17_trunc_ratio, group=="missense in ND4"), aes(x=gene,y=log_FC), size=3.5, pch=18) + 
  scale_color_cosmic(palette="signature_substitutions") + 
  theme_classic() + coord_cartesian(ylim=c(-0.75,1.25))




### Fig6e. DNA & RNA VAF compare -----------------------------------

rRNA_struc1 <- read_tsv("~/reference/MT/mtRNA/d.16.m.H.sapiens.5.bpseq_.txt", col_names = "info") %>% 
  mutate(region2="MT-RNR1")
rRNA_struc2 <- read_tsv("~/reference/MT/mtRNA/d.235_233.m.H.sapiens1.bpseq.txt", col_names="info") %>% 
  mutate(region2="MT-RNR2")
rRNA_struc <- rbind(rRNA_struc1,rRNA_struc2) %>% 
  separate(info,into=c("rnaPOS","REF","pairing"), sep=" ") %>% 
  mutate(across(c(rnaPOS,pairing), as.numeric))

rRNA_struc_mod <- rRNA_struc %>% 
  mutate(start=ifelse(region2=="MT-RNR1",647,1670), end=ifelse(region2=="MT-RNR1",1601,3229)) %>% 
  mutate(POS=start+rnaPOS)  %>% 
  mutate(POS=ifelse(POS>=3107,POS+1,POS)) %>% 
  mutate(rna_region=ifelse(pairing>0,"stem","loop"))

tRNA_struc <- read_tsv("~/reference/MT/refGene_MTs_tRNAstructure.txt")
tRNA_struc_mod <- tRNA_struc %>% 
  select(gene,gstart,gstop,ACCstem1:Terminal) %>% 
  gather(4:ncol(.),key="region",value="len") %>% 
  arrange(gene) %>% 
  group_by(gene) %>% 
  mutate(cumsum=cumsum(len)) %>% 
  mutate(nextpos=lag(cumsum)) %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  mutate(start=gstart+nextpos+1,end=gstart+cumsum) %>% 
  rowwise() %>% 
  mutate(POS=paste(seq(start,end),collapse = ":")) %>% 
  separate_rows(POS,sep=":") %>% 
  mutate(rna_region=ifelse(grepl("stem",region),"stem",ifelse(grepl("link",region),"link", "loop")))

MT_rna_struc <- rRNA_struc_mod %>% 
  dplyr::rename("gene"=region2) %>% 
  select(gene,POS,rna_region) %>% 
  rbind(tRNA_struc_mod %>% select(gene,POS,rna_region))


# not tRNA & rRNA region
df_rna_vaf %>% 
  filter(TD>100) %>% 
  filter(consequence %in% c("intergenic","synonymous","missense","truncating")) %>% 
  filter(region!=".") %>% 
  mutate(group=ifelse(consequence=="intergenic",region, consequence)) %>% 
  ggplot() + 
  geom_point(aes(x=VAF,y=rnaVAF)) + 
  geom_abline(intercept=1, col="gray20") + 
  facet_wrap(~factor(consequence,levels=c("intergenic","synonymous","missense","truncating")), scales="free") + 
  coord_cartesian(xlim=c(0,100),ylim=c(0,100)) +
  theme_classic() + 
  theme(strip.text.x = element_text(color="black")) 

# tRNA
df_rna_vaf %>% filter(vartype!="fe") %>% filter(vartype_2!="frequent_recurrent") %>% 
  filter(TD>100) %>% 
  filter(region=="tRNA") %>% 
  merge(tRNA_struc_mod,by=c("POS")) %>% 
  as_tibble() %>% 
  ggplot() + 
  geom_abline(intercept=1, col="gray20")+
  geom_point(aes(x=VAF,y=rnaVAF,col=rna_region)) + 
  scale_color_jama() + 
  theme_classic() +
  coord_cartesian(xlim=c(0,100),ylim=c(0,100)) 


# rRNA
df_rna_vaf %>% filter(vartype!="fe") %>% filter(vartype_2!="frequent_recurrent") %>% 
  filter(TD>100) %>% 
  filter(region=="rRNA") %>% 
  merge(rRNA_struc_mod,by=c("POS")) %>% 
  as_tibble() %>% 
  ggplot() + 
  geom_abline(intercept=1, col="gray20")+ 
  geom_point(aes(x=VAF,y=rnaVAF, col=rna_region)) + 
  theme_classic() + 
  theme(strip.text.x = element_text(color="black")) + 
  scale_color_jama()+
  coord_cartesian(xlim=c(0,100),ylim=c(0,100))




### Fig5f. tRNA region ---------------------------------------------

# tRNA
df_rna_vaf_trna <- df_rna_vaf %>% filter(vartype!="fe") %>% filter(vartype_2!="frequent_recurrent") %>% 
  filter(TD>100) %>% 
  filter(region=="tRNA") %>% 
  mutate(diff=rnaVAF/VAF) %>% 
  merge(MT_rna_struc,by=c("POS")) %>% 
  as_tibble() %>% 
  filter(rna_region!="link") 

df_rna_vaf_trna %>% 
  ggplot(aes(x=rna_region,y=log2(diff))) + 
  geom_jitter(aes(col=rna_region),width=0.1) + 
  geom_violin(alpha=0)  + 
  scale_color_jama() + theme_classic() + ggtitle("tRNA")

wilcox.test(df_rna_vaf_trna %>% filter(rna_region=="stem") %>% pull(diff), 
  df_rna_vaf_trna %>% filter(rna_region=="loop") %>% pull(diff), alternative = "greater")  #pval=0.01578



### Fig5g. rRNA region ---------------------------------------------

# rRNA
df_rna_vaf_rrna <- df_rna_vaf %>% filter(vartype!="fe") %>% filter(vartype_2!="frequent_recurrent") %>% 
  filter(TD>100) %>%
  filter(region=="rRNA") %>% 
  mutate(diff=rnaVAF/VAF) %>% 
  mutate(diff2=rnaVAF-VAF) %>% 
  merge(MT_rna_struc,by=c("POS")) %>% 
  as_tibble() %>% 
  filter(rna_region!="link") 

df_rna_vaf_rrna %>% 
  ggplot(aes(x=rna_region,y=log2(diff))) + 
  geom_jitter(aes(col=rna_region),width=0.1) + 
  geom_violin(alpha=0)  + 
  scale_color_jama() + theme_classic() + ggtitle("rRNA")

wilcox.test(df_rna_vaf_rrna %>% filter(rna_region=="stem") %>% pull(diff), 
  df_rna_vaf_rrna %>% filter(rna_region=="loop") %>% pull(diff), alternative = "less")  #pval=0.0329