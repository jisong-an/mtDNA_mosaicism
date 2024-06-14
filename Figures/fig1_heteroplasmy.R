################################
#
#   Figure 1 -mtDNA heteroplasmy in normal clones
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


### Fig1b. mutational landscape -------------------------------------------------------------------

plot_vargraph <- function(temp_db, max, tree_list){
  
  tot_pat <- nrow(temp_db)
  pat_pos <- temp_db %>% select(patient)  %>% select(patient) %>% group_by(patient) %>% 
    mutate(n=n()) %>% ungroup() %>% unique() %>% pull(n)
  pat_name <- temp_db %>% select(patient) %>% unique() %>% pull(patient)
  pat_age <- temp_db %>% select(patient,age) %>% unique() %>% pull(age)
  
  # x axis -patient
  xaxis_pos = c(0)
  for(i in 1:length(pat_pos)){xaxis_pos[i+1] = sum(pat_pos[1:i])}
  xaxis_plot_pos = xaxis_pos/tot_pat
  
  grid.newpage()
  
  ## 1. Phylogeny  [0.85-1.00]
  pushViewport(viewport(x=0, y=0, width=1, height=1, just=c('left','bottom')))
  pushViewport(viewport(x=0, y=0.85, height=0.1, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height=1, width=0.9))
  
  for (i in 1:length(pat_pos)){
    pushViewport(viewport(x=xaxis_plot_pos[i], y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1]-xaxis_plot_pos[i]))  
    pat = pat_name[i]
    proj = temp_db %>% filter(patient==pat) %>% pull(project) %>% .[1]
    if (proj=="DB") {f=grep(paste0("DB",pat), tree_list, value=TRUE)}
    else {f=grep(pat,tree_list,value=TRUE)}
    tree <- read.tree(f)
    gt <- ggtree(tree, branch.length="none", size=0.1) + layout_dendrogram() + 
      theme(plot.margin = unit(c(0,0,0,0),"cm"))
    g <- ggplotGrob(gt)
    grid.draw(g)
    popViewport(1)
  }
  
  popViewport(3)
  
  
  ## 2. Total variant count [0.40-0.85]
  pushViewport(viewport(x=0, y=0.4, height=0.45, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height = 1, width=0.9))
  grid.rect(gp=gpar(lty="solid"))
  
  data_n <- temp_db
  data_order <- data_n %>% mutate(order=row_number()/tot_pat) %>% pull(order)
  
  # y-value : variant count
  yrange=max
  yaxis_plot_pos = seq(0,1,by=1/yrange*4)
  y_value <- data_n %>% mutate(n=n/yrange) %>% pull(n)
  y_value2 <- data_n %>% mutate(feN=feN/yrange) %>% pull(feN)
  
  # x : each clone
  xaxis_data <- data_n %>% mutate(order=row_number()/tot_pat) %>% pull(order)
  
  # background
  for (i in 1:length(pat_pos)){
    if(i %% 2 == 0){bgcolor="gray92"}
    else{bgcolor = FALSE}
    grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
  }
  
  # axis
  grid.yaxis(at=yaxis_plot_pos, label=yaxis_plot_pos*yrange, gp = gpar(fontsize=7))
  grid.text('total variant count', x=-0.04, rot=90, gp=gpar(fontsize=9))
  
  
  # bar graph
  for (i in 1:tot_pat){
    grid.rect(x=xaxis_data[i],y=0, just=c('right','bottom'), height=y_value[i], width=xaxis_data[1], gp=gpar(fill="#5B6777", col=FALSE))  # total variant
    grid.rect(x=xaxis_data[i],y=0, just=c('right','bottom'), height=y_value2[i], width=xaxis_data[1], gp=gpar(fill="#A0C7EA", col=FALSE))  # fe variant
  }
  
  popViewport(2)
  
  
  ## 3. maximum VAF [0.10-0.40]
  pushViewport(viewport(x=0, y=0.1, height=0.3, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height = 0.9, width=0.9))
  grid.rect(gp=gpar(lty="solid"))
  
  # background
  for (i in 1:length(pat_pos)){
    if(i %% 2 == 0){bgcolor="gray92"}
    else{bgcolor = FALSE}
    grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
  }
  
  # axis
  grid.yaxis(at=c(0,0.25,0.5,0.75,1), label=c(0,25,50,75,100), gp = gpar(fontsize=7))
  grid.text('maximum VAF (%) ', x=-0.04, rot=90, gp=gpar(fontsize=9))
  
  # y-value : max VAF
  y_value_vaf <- data_n %>% mutate(maxVAF = maxVAF/100) %>% pull(maxVAF)
  
  # bar graph
  for (i in 1:tot_pat){
    grid.rect(x=xaxis_data[i],y=0, just=c('right','bottom'), height=y_value_vaf[i], width=xaxis_data[1], gp=gpar(fill="lightcyan4", col=FALSE))
    #grid.rect(x=xaxis_data[i],y=1, just=c('right','top'), height=1-y_value_vaf[i], width=xaxis_data[1], gp=gpar(fill="gray", col=FALSE))
  }
  
  popViewport(2)
  
  
  ## 4. tissue [0.05-0.1]
  pushViewport(viewport(x=0,y=0.05, height=0.05, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))
  
  grid.text('tissue', x=-0.03, gp=gpar(fontsize=8))
  
  for (i in 1:length(pat_pos)){
    pat = pat_name[i]
    proj = temp_db %>% filter(patient==pat) %>% pull(project) %>% .[1]
    if (proj=="DB") {bgcolor=pal_mt[[2]]}
    else if (proj=="Line1" | proj=="Abortus") { bgcolor=pal_mt[[1]]}
    else{bgcolor=pal_mt[[3]]}
    grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=bgcolor, col=FALSE))
    #grid.text(pat_name[i], x=(xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=-0.04, gp=gpar(fontsize=7.5))
  }
  popViewport(2)
  
  
  ## 5. age [0-0.05]
  pushViewport(viewport(x=0,y=0, height=0.05, just=c('left','bottom')))
  pushViewport(viewport(x=0.5, y=0.5, height=0.9, width=0.9))
  
  grid.text('age', x=-0.03, gp=gpar(fontsize=8))
  col_fun = colorRamp2(c(0,30,100),c("lightgoldenrod","goldenrod1","darkorchid4"))
  #col_fun = colorRamp2(c(0,30,100),c("#FDFBF5","#FF97B5","#6E2A85"))
  
  for (i in 1:length(pat_pos)){
    grid.rect(x=xaxis_plot_pos[i],y=0, just=c('left','bottom'), height=1, width=xaxis_plot_pos[i+1] - xaxis_plot_pos[i], gp=gpar(fill=col_fun(pat_age[i]), col=FALSE))
    #grid.text(pat_name[i], x=(xaxis_plot_pos[i]+xaxis_plot_pos[i+1])/2, y=-0.04, gp=gpar(fontsize=7.5))
  }
  popViewport(2)
  
}

tree_file <- c("~/project/11_Clone_MT/info/11_Abortus/tree/SA03.length.nwk", getAbsolutePath("~/project/11_Clone_MT/12_Heatmap/01_Line1/","*nwk"), getAbsolutePath("~/project/11_Clone_MT/12_Heatmap/02_DB/","*nwk"), "~/project/11_Clone_MT/info/04_H.Blood/tree/BM_tree.nwk",getAbsolutePath("~/project/11_Clone_MT/info/09_blood_new/tree/","*nwk"))
data_sample_order <- tibble()

for (f in tree_file){
  pat <- gsub("DB","",strsplit(basename(f), split="\\.|_")[[1]][1])
  tree <- read.tree(f)
  print(ggtree(tree, branch.length="none") + layout_dendrogram())
  sample_order <- fortify(tree) %>% filter(!is.na(label)) %>% arrange(y) %>% filter(isTip==TRUE) %>% pull(label)
  if (pat %in% c("10","2","3","5","6","8","9")){
    meta_dt <- read_tsv("/home/users/anjisong/project/11_Clone_MT/12_Heatmap/02_DB/Summary_per_sample_210330.txt") %>% 
      filter(deadbody==paste0("DB",pat) & current_final_for_lineage == 'Y') %>% 
      select(lineage_id, sample_id)
    sample_order <- fortify(tree) %>% merge(meta_dt %>% rename("label"=lineage_id), by=c("label"), all.x = TRUE) %>% 
      as_tibble() %>% filter(!is.na(sample_id)) %>% arrange(y) %>% pull(sample_id)
  }
  data_sample_order <- rbind(data_sample_order,sample_info_final_v3 %>% filter(patient==pat) %>% filter(type!="bulk_normal") %>% arrange(factor(sample,levels=sample_order)) %>% select(patient,sample))  
}

temp <- df %>% 
  filter(!is.na(VAR)) %>% filter(type=="clonal_normal") %>% 
  group_by(project, patient, VAR) %>% 
  mutate(varn=n()) %>% select(varn,everything()) %>% 
  mutate(vartype_3=ifelse(varn==1,"unique","shared")) %>% 
  group_by(project,patient,sample) %>% 
  mutate(feN=sum(vartype_3=="shared")) %>% 
  mutate(maxVAF=max(VAF)) %>% 
  group_by(project,patient,sample,tissue2, age, maxVAF,feN) %>% 
  tally() %>% 
  merge(sample_info_final_v3 %>% filter(type=="clonal_normal") %>% select(project,patient,sample,age), by=c("project","patient","sample","age"), all=TRUE) %>% 
  as_tibble() %>% 
  filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  mutate_all(~replace(., is.na(.),0)) %>% arrange(factor(sample, levels=data_sample_order %>% pull(sample)))


pdf("~/project/11_Clone_MT/figure/variant_maxVAF_dist.240216_shared_tree.pdf", width=20,height=4.5, onefile = T)
p1 <- plot_vargraph(temp %>% filter(project %in% c("Abortus","Line1","DB", "Hblood")) %>% 
                      arrange(factor(project, levels=c("Abortus","Line1","DB","Hblood")),age), 
                    max=15, tree_list=tree_file)
p2 <- plot_vargraph(temp %>% filter(project %in% c("blood_new")) %>% arrange(age), max=15, tree_list=tree_file)
dev.off()



### Fig1c. Circos --------------------------------------------------

draw_MT_circos_scale <- function(df){
  
  mito_genes <- as_tibble(read.table("/home/users/kimin/projects/12_MT/00_reference/mito_genes.bed", col.names = c("chrom","start","end","gene","dummy","strand","type"), colClasses = c("character","integer","integer","character","integer","character","character")))
  mito_genes$direction <- ifelse(mito_genes$strand == "+", 1, -1)
  mito_genes$angle <- 360 * (mito_genes$start + mito_genes$end)/2/16569
  mito_genes$gene_label <- ifelse((mito_genes$end-mito_genes$start) < 350, "", mito_genes$gene)
  
  #gene_color_vector <- brewer.pal(7, "Set3")
  gene_color_vector <- pal_jama()(7)
  names(gene_color_vector) <- sort(unique(mito_genes$type))
  
  snp_color_vector = c("C>A"="#15A0EC","G>T"="#15A0EC", # blue
                       "C>G"="#0D0C1B","G>C"="#0D0C1B", # black
                       "C>T"="#F23A29","G>A"="#F23A29", # red
                       "T>A"="#A1A1A1","A>T"="#A1A1A1", # grey
                       "T>C"="#5AB440","A>G"="#5AB440", # green
                       "T>G"="#F2BBC5","A>C"="#F2BBC5") # pink
  
  snp_pch_vector = c("intergenic"=15,                 # square
                     "synonymous"=16,                 # circle
                     "inframe indel"=18,              # diamond
                     "missense"=17,                   # triangle
                     "truncating"=4,                  # x
                     "frameshift"=6,                  # empty triangle
                     "incomplete_terminal_codon"=7)   # empty square x
  
  genes <- ggplot(mito_genes, aes(xmin=start, xmax=end, y=chrom, fill=type, label=gene, forward=direction)) + 
    geom_gene_arrow(arrowhead_height=unit(3, "mm"), arrowhead_width=unit(1, "mm")) +
    geom_gene_label(align = 'centre') +
    facet_wrap(~ chrom, scales='free', ncol=1) +
    scale_fill_jama() +
    theme_genes() + guides(fill="none") + 
    theme(legend.position='bottom') + ylab('') + xlab('') + 
    scale_x_continuous(limits = c(0,16569), expand = c(0, 0), breaks = seq(0, 16000, by = 4000))
  
  
  
  circos.clear()
  circos.par(start.degree=90, gap.degree=1, cell.padding=c(0.01,0,0.01,0), xaxis.clock.wise=FALSE)
  circos.initialize(sectors=factor(c("MT")), xlim=c(0,16569))
  
  
  # Heavy strand
  
  #df_heavy <- df %>% filter(strand=="H") %>% mutate(VAFln=ifelse(VAF<5, VAF*10, (VAF-5)/95*50+50))
  df_heavy <- df %>% filter(strand=="H") %>% mutate(VAFln=ifelse(VAF<10, VAF*5, (VAF-10)/90*50+50))
  
  circos.track(ylim=c(0,100), track.height=0.2, cell.padding=c(0,0,0,0), bg.border=NA)
  circos.lines(c(0,16569), c(100,100), col='gray')
  circos.lines(c(0,16569), c(50,50), col='gray')
  circos.lines(c(0,16569), c(0,0), col='gray')
  
  circos.points(x=df_heavy$POS, y=df_heavy$VAFln, col=snp_color_vector[df_heavy$var],
                pch=snp_pch_vector[df_heavy$consequence], cex = 0.8)
  circos.yaxis(side='right',at=c(0,50,100), labels=c("0","0.1","1"), tick=TRUE, labels.cex = 0.7, tick.length = 0.3)
  set_track_gap(0.06)
  
  
  
  
  # MITOMAP
  circos.track(ylim=c(0,1),track.height=0.1, bg.border=NA)
  for(i in 1:nrow(mito_genes)){
    
    if( mito_genes[[i,'strand']] == '+'){
      
      circos.arrow(x1=mito_genes[[i,'start']], x2=mito_genes[[i,'end']], y=0.5, width=1, 
                   arrow.head.width = 1, arrow.head.length = 50, arrow.position = 'end', 
                   col=gene_color_vector[mito_genes[[i,'type']]], border="black")
      
    } else if (mito_genes[[i,'strand']] == '-'){
      
      circos.arrow(x1=mito_genes[[i,'start']], x2=mito_genes[[i,'end']], y=0.5, width=1, 
                   arrow.head.width = 1, arrow.head.length = 50, arrow.position = 'start', 
                   col=gene_color_vector[mito_genes[[i,'type']]], border="black")
      
    } else {
      print("ERROR")
    }
    
  }
  
  circos.text(x=(mito_genes$start+mito_genes$end)/2, y=rep(0.5, nrow(mito_genes)), labels=mito_genes$gene_label, facing='inside',cex = 0.8)
  circos.xaxis(h='top',major.at=c(seq(0,16569,2000), 16569), labels=c(seq(0,16569,2000)), minor.ticks=3, labels.cex=0.8, major.tick.length=0.1)
  
  set_track_gap(0.03)
  
  
  # Light strand
  
  #df_light <- df %>% filter(strand=="L") %>% mutate(VAFln=ifelse(VAF<5, VAF*10, (VAF-5)/95*50+50))
  df_light <- df %>% filter(strand=="L") %>% mutate(VAFln=ifelse(VAF<10, VAF*5, (VAF-10)/90*50+50))
  
  circos.track(ylim=c(0,100), track.height=0.2, cell.padding=c(0,0,0,0), bg.border=NA)
  circos.lines(c(0,16569), c(100,100), col='gray')
  circos.lines(c(0,16569), c(50,50), col='gray')
  circos.lines(c(0,16569), c(0,0), col='gray')
  
  circos.points(x=df_light$POS, y=df_light$VAFln, col=snp_color_vector[df_light$var],
                pch=snp_pch_vector[df_light$consequence], cex=0.8)
  
  circos.yaxis(side='right',at=c(0,50,100), labels=c("0","0.1","1"), tick=TRUE, labels.cex = 0.7, tick.length = 0.3)
  
  
}


df_sofe <- df  %>% filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(type=="clonal_normal") %>% 
  filter(!is.na(VAR)) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(var=paste0(REF,">",ALT)) %>%
  mutate(strand=ifelse(REF=="C" | REF=="T", "L", "H")) %>% 
  mutate(var=ifelse(var=="G>A","C>T",ifelse(var=="G>T","C>A", ifelse(var=="G>C","C>G",ifelse(var=="A>T","T>A",ifelse(var=="A>G","T>C",ifelse(var=="A>C","T>G",var))))))) %>% 
  filter(consequence %in% c("intergenic","missense","synonymous","truncating")) %>% 
  mutate(vartype=ifelse(vartype=="fe","fe",ifelse(vartype %in% c("","lineage","recurrent"),"somatic",vartype))) %>%
  mutate(VAF2=ifelse(vartype=="fe",pseudoVAF,VAF)) %>% 
  select(-VAF,-sample, -mtDNA_depth,-nDNA_depth,-mtCN) %>% 
  unique() %>% rename("VAF"=VAF2)

draw_MT_circos_scale(df_sofe)



### Fig1d. Signature -------------------------------------------------

df %>% 
  filter(type=="clonal_normal") %>% 
  filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>%  
  select(project,patient,VAR,POS:Gene,tissue2) %>% 
  unique() %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(var=paste0(REF,">",ALT)) %>%
  mutate(strand=ifelse(REF=="C" | REF=="T", "L", "H")) %>% 
  mutate(sigvar=ifelse(var=="G>A","C>T",ifelse(var=="G>T","C>A", ifelse(var=="G>C","C>G",ifelse(var=="A>T","T>A",ifelse(var=="A>G","T>C",ifelse(var=="A>C","T>G",var))))))) %>% 
  group_by(tissue2,sigvar,strand) %>% tally() %>%
  group_by(tissue2) %>%  
  mutate(percent=n/sum(n)) %>% 
  ggplot(aes(x=sigvar,y=percent, fill=sigvar, alpha=strand)) + 
  geom_hline(yintercept=c(0,0.2,0.4),color="gray90") + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~factor(tissue2,levels=c("colon","fibroblast","blood"))) + 
  theme_classic() +
  #geom_text(aes(label=round(percent,1)), position=position_dodge(width=1), vjust=-0.5) + 
  ggtitle("mutational signature") + 
  scale_fill_cosmic(palette = "signature_substitutions") + 
  scale_alpha_manual(values=c(1,0.5)) 



### Fig1e. strand-bias ------------------------------------------

## Count
df_sig <- df  %>% filter(project %in% c("Line1","DB","Hblood","blood_new","Abortus")) %>% 
  filter(type=="clonal_normal") %>% 
  filter(!is.na(VAR)) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(var=paste0(REF,">",ALT)) %>%
  mutate(strand=ifelse(REF=="C" | REF=="T", "L", "H")) %>% 
  mutate(var=ifelse(var=="G>A","C>T",
                    ifelse(var=="G>T","C>A", 
                           ifelse(var=="G>C","C>G",
                                  ifelse(var=="A>T","T>A",
                                         ifelse(var=="A>G","T>C",
                                                ifelse(var=="A>C","T>G",var)))))))

df_sig %>% 
  mutate(regionRep=ifelse(POS<192 | POS > 16196, "RepOrigin","OutRep")) %>% 
  filter(var %in% c("C>T","T>C")) %>% 
  group_by(regionRep,var,strand) %>% 
  tally() %>% mutate(REF=str_remove(var,">.*")) %>% 
  merge(reference %>% 
          select(1:3) %>% 
          unique() %>% 
          mutate(regionRep=ifelse(POS<192 | POS > 16196, "RepOrigin","OutRep")) %>% 
          group_by(REF,regionRep) %>% 
          tally() %>% 
          mutate(strand=ifelse(REF %in% c("C","T"),"L","H")) %>% 
          mutate(REF=ifelse(REF=="A","T", ifelse(REF=="G","C",REF))) %>% 
          dplyr::rename("background"=n),by=c("REF","regionRep","strand"),all=TRUE) %>% 
  as_tibble() %>% 
  mutate(normalize_n=n/background) %>% 
  select(REF,regionRep,var,strand,normalize_n) %>% 
  spread(key="strand",value="normalize_n") %>% mutate(bias=H/L) %>% 
  ggplot() + 
  geom_bar(aes(x=var,y=log2(bias),fill=regionRep), position=position_dodge(), stat="identity") + 
  xlab("") + 
  theme_classic() + 
  scale_fill_brewer(palette = "Accent")



### Fig1f. hypermutated ----------------------
  
  mito_genes <- as_tibble(read.table("/home/users/kimin/projects/12_MT/00_reference/mito_genes.bed", col.names = c("chrom","start","end","gene","dummy","strand","type"), colClasses = c("character","integer","integer","character","integer","character","character")))
  mito_genes$direction <- ifelse(mito_genes$strand == "+", 1, -1)
  mito_genes$gene_label <- ifelse((mito_genes$end-mito_genes$start) < 350, "", mito_genes$gene)
  
  gene_color_vector <- pal_jama()(7)
  names(gene_color_vector) <- sort(unique(mito_genes$type))
  
  genes <- ggplot(mito_genes, aes(xmin=start, xmax=end, y=chrom, fill=type, label=gene, forward=direction)) + 
    geom_gene_arrow(arrowhead_height=unit(3, "mm"), arrowhead_width=unit(1, "mm")) +
    geom_gene_label(align = 'centre') +
    facet_wrap(~ chrom, scales='free', ncol=1) +
    scale_fill_jama()+
    theme_genes() + guides(fill="none") + 
    theme(legend.position='bottom') + ylab('') + xlab('') + 
    scale_x_continuous(limits = c(0,16569), expand = c(0, 0), breaks = seq(0, 16000, by = 4000))
  
  
  ## 10_ARL10-4_001D4
  g <- df %>% filter(sample=="10_ARL10-4_001D4") %>% mutate(var=paste0(REF,">",ALT)) %>%
    mutate(strand=ifelse(REF=="C" | REF=="T", "L", "H")) %>% 
    mutate(sigvar=ifelse(var=="G>A","C>T",ifelse(var=="G>T","C>A", ifelse(var=="G>C","C>G",ifelse(var=="A>T","T>A",ifelse(var=="A>G","T>C",ifelse(var=="A>C","T>G",var))))))) %>% 
    ggplot(aes(x=POS,y=VAF,col=sigvar)) + 
    geom_jitter(height=0.1, size=2.5) + 
    scale_y_continuous(limits=c(0,15)) + 
    scale_x_continuous(limits = c(0,16569), expand=c(0,0), breaks=seq(0,16000,by=4000)) + 
    scale_color_manual(values = pal_cosmic("signature_substitutions")(6)[c(3,4,5,6)]) + 
    theme_classic() + 
    ylab("VAF (%)") + 
    ggtitle("10_ARL10-4_001D4 (normal fibroblast) ") + 
    theme(legend.position=c(0.95,0.7)) 
  
  
  plot_grid(g,NULL,genes,ncol=1, rel_heights=c(1.5,-0.02,0.5))
