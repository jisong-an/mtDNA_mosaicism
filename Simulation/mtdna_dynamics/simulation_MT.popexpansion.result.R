################################
#
#   mtDNA expansion simulation result
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################


library(tidyverse)
library(RColorBrewer)
library(circlize)


### All simulations -----------------------------------------------

input <- getAbsolutePath("~/project/11_Clone_MT/10_Simulation/231229_popexpansion/","*tsv")

pdf("~/project/11_Clone_MT/10_Simulation/simul_population_expansion.240110.v2.pdf", width=10, height=6, onefile=T)
for (f in input){
  
  cn <- gsub("mtcn","",strsplit(basename(f),split = "\\.")[[1]][2])
  temp <- read_tsv(f) %>% dplyr::filter(generation <= 3000)
  temp <- temp %>% 
    gather(starts_with("cell"), key="cell",value="VAFlist") %>% 
    separate(VAFlist, sep="\\|", into=c("VAFlist","x")) %>% 
    select(-x) %>% 
    mutate(VAFlist=ifelse(VAFlist=="",0,VAFlist)) %>% 
    separate(VAFlist, sep=":", into=c(paste0("VAF",seq(1:10)))) %>% 
    mutate_all(~replace(.,is.na(.),0))  %>% 
    mutate(across(VAF1:VAF10,as.numeric)) %>% 
    mutate(sumVAF=VAF1+VAF2+VAF3+VAF4+VAF5+VAF6+VAF7+VAF8+VAF9+VAF10) %>% 
    mutate(nVAR=rowSums(across(VAF1:VAF10)>0))
  temp_mod <- temp %>% dplyr::select(generation,mtCN,cell,VAF1) %>% spread(key="cell",value="VAF1")
  
  assign(paste0("mtDNA_pop_expansion_cn",cn), temp)
  assign(paste0("mtDNA_pop_expansion_cn",cn,"_mod"), temp_mod)
  
  print(temp %>% 
          dplyr::filter(generation<=1500) %>% 
          ggplot() + 
          stat_density_2d(geom="polygon", 
                          mapping=aes(x=generation, y=VAF1, fill=after_stat(level)),
                          contour=T, size=0.15,contour_var = "ndensity", 
                          color="white", bins=20)+
          scale_x_continuous(expand=c(0,0), limits=c(0,1550), breaks=seq(0,1500,300)) +
          scale_y_continuous(expand=c(0,0), limits=c(0,1.1)) + 
          scale_fill_gradientn(colors=c("#fdf0f0","#d6c9c9","#c5aec9","#9e8fa0","#918493","#8079a9" ,"#78709b","#656083","#535172","#302c3e"))+
          theme_classic() +
          theme(text=element_text(size=15), legend.title=element_blank(), 
                panel.border = element_rect(color="black", fill=NA)) + 
          xlab("turnover") + ylab("VAF") + ggtitle(paste0("mtDNA copy number = ",cn))
  )  
}

dev.off()

save(list=c(ls()[grepl("^mtDNA_pop", ls())]), file="~/R/R_script/11_Clone_MT/mtDNA_population_expansion_simul.all.240110.RData")