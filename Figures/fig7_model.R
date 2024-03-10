################################
#
#   Figure 7 -mtDNA dynamics
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

library(tidyverse)
library(RColorBrewer)
library(grid)
library(circlize)


### Fig7b. mtDNA expansion (mtCN=750) --------------------------------------

  load("~/R/R_script/11_Clone_MT/mtDNA_population_expansion_simul.all.240110.RData")

  mtDNA_pop_expansion_cn750 %>% 
    dplyr::filter(generation<=1500) %>% 
    filter(generation %% 10 == 0) %>% 
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
    xlab("Mitotic turnover") + ylab("VAF") + ggtitle("mtDNA copy number = 750")



### Fig7c. mtDNA heteroplasmy expansion (mtCN=750) --------------------------

  load("~/R/R_script/11_Clone_MT/mtDNA_mutation_expansion_simul.all.mod2.240131.RData")
  
  mtDNA_mut_expansion_cn750 %>% 
    dplyr::filter(generation<=1500) %>% 
    filter(generation %% 10 == 0) %>% 
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
    xlab("Mitotic turnover") + ylab("VAF") + ggtitle("mtDNA copy number = 750")