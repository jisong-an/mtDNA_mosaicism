#######################
#
#  Simulation: HetFE variant cavaf sampling
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
#######################

# 1. sampling cells from the simulation results (simulation_MT.fe_cavaf.<turnover>.py) to mimic sequencing
# 2. calculate summary statistics in observed data
# 3. calculate MSE (mean squared error) between simulated & observed data

# modification of simulation directory & Rdata name is needed


suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(magrittr))
suppressMessages(library(tidyverse))


## Data
load("/home/users/anjisong/R/R_script/11_Clone_MT/231111_Rdata_workspace.RData")


## Functions
getAbsolutePath <- function(dir, pattern) {
    
    return(paste0(normalizePath(dir), "/", list.files(normalizePath(dir), pattern)))
    
  }

filter <- dplyr::filter
select <- dplyr::select

pharse_obs <- function(df,pat,var,remove_sample){
  
  vafdist_col <- c(seq(5,100,by=5) %>% as.character())
  vafrange_df <- as_tibble(matrix(nrow = 0, ncol = length(vafdist_col)), 
                           .name_repair = ~ vafdist_col)
  obs_temp <- df %>% 
    filter(VAR==var) %>% 
    filter(patient==pat) %>% 
    dplyr::select(patient,sample,VAR,VAF,age,tissue2) %>% 
    merge(sample_info_final_v3 %>% filter(type=="clonal_normal") %>% 
            filter(patient==pat) %>% filter(!(sample %in% remove_sample)) %>% 
            dplyr::select(patient,sample,age,tissue2), by=c("patient","sample","age","tissue2"), all=TRUE) %>% 
    as_tibble() %>% 
    mutate(across(VAR,~replace(.,is.na(.),var)), across(VAF,~replace(.,is.na(.),0))) 
  
  obs_temp_mod <- obs_temp %>% 
    dplyr::select(sample,VAF) %>% 
    mutate(wt_ratio=sum(ifelse(VAF<0.5, 1, 0))/nrow(obs_temp), 
           hetero_ratio=sum(ifelse(VAF>=0.5 & VAF<90, 1, 0))/nrow(obs_temp), 
           homo_ratio=sum(ifelse(VAF>=90, 1, 0))/nrow(obs_temp)) %>% 
    mutate(VAF=VAF/100) %>% 
    mutate(meanVAF=round(mean(VAF),3), sdVAF=round(sd(VAF),3)) %>% 
    mutate(x=((VAF*100)%/%5)+1) %>% 
    mutate(group=ifelse(VAF==0,0,x*5)) %>% 
    group_by(group,wt_ratio,hetero_ratio,homo_ratio,meanVAF,sdVAF) %>% 
    tally() %>% 
    filter(group!=0) %>% 
    mutate(n=n/nrow(obs_temp)) %>% 
    spread(key="group",value="n") %>% 
    merge(vafrange_df, all=TRUE) %>% 
    as_tibble() %>% 
    mutate_all(~replace(.,is.na(.),0)) %>% 
    dplyr::select(wt_ratio:sdVAF,vafdist_col)
  
  invisible(obs_temp_mod)
  
}

calculate_distance <- function(simul,obs){
  
  simul_res <- simul %>%
    mutate(distance=((wt_ratio-obs[["wt_ratio"]])^2 + 
                       (hetero_ratio-obs[["hetero_ratio"]])^2 + 
                       (homo_ratio-obs[["homo_ratio"]])^2 + 
                       (meanVAF-obs[["meanVAF"]])^2 + 
                       (sdVAF-obs[["sdVAF"]])^2 + 
                       (`5`-obs[["5"]])^2+ 
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
                       (`100`-obs[["100"]])^2)) 
  
  invisible(simul_res)
  
}
  

sample_cells_sumstat <- function(df,start,end,patN, iter){
    simul_merge <- tibble()
    for (i in 1:iter){
      print(paste0(i," iteration start!"))
      idx = sample(seq(start,end),size=patN,replace=FALSE)
      temp <- df %>% 
        select(idx,generation) %>% 
        mutate(wt_ratio=rowSums(across(1:patN)<0.005)/patN, 
               hetero_ratio=rowSums(across(1:patN)>=0.005 & across(1:patN)<0.90)/patN, 
               homo_ratio=rowSums(across(1:patN)>=0.90)/patN) %>% 
        mutate(`5`=rowSums(across(1:patN)>=0.005 & across(1:patN)<0.05)/patN, 
               `10`=rowSums(across(1:patN)>=0.05 & across(1:patN)<0.1)/patN,
               `15`=rowSums(across(1:patN)>=0.1 & across(1:patN)<0.15)/patN,
               `20`=rowSums(across(1:patN)>=0.15 & across(1:patN)<0.2)/patN,
               `25`=rowSums(across(1:patN)>=0.2 & across(1:patN)<0.25)/patN,
               `30`=rowSums(across(1:patN)>=0.25 & across(1:patN)<0.3)/patN,
               `35`=rowSums(across(1:patN)>=0.3 & across(1:patN)<0.35)/patN,
               `40`=rowSums(across(1:patN)>=0.35 & across(1:patN)<0.4)/patN,
               `45`=rowSums(across(1:patN)>=0.4 & across(1:patN)<0.45)/patN,
               `50`=rowSums(across(1:patN)>=0.45 & across(1:patN)<0.5)/patN,
               `55`=rowSums(across(1:patN)>=0.5 & across(1:patN)<0.55)/patN,
               `60`=rowSums(across(1:patN)>=0.55 & across(1:patN)<0.6)/patN,
               `65`=rowSums(across(1:patN)>=0.6 & across(1:patN)<0.65)/patN,
               `70`=rowSums(across(1:patN)>=0.65 & across(1:patN)<0.7)/patN,
               `75`=rowSums(across(1:patN)>=0.7 & across(1:patN)<0.75)/patN,
               `80`=rowSums(across(1:patN)>=0.75 & across(1:patN)<0.8)/patN,
               `85`=rowSums(across(1:patN)>=0.8 & across(1:patN)<0.85)/patN,
               `90`=rowSums(across(1:patN)>=0.85 & across(1:patN)<0.9)/patN,
               `95`=rowSums(across(1:patN)>=0.9 & across(1:patN)<0.95)/patN,
               `100`=rowSums(across(1:patN)>=0.95 & across(1:patN)<1.1)/patN) %>% 
        rowwise() %>% 
        mutate(meanVAF=mean(c_across(starts_with("cell"))), sdVAF=sd(c_across(starts_with("cell")))) %>% 
        select(generation,wt_ratio,hetero_ratio,homo_ratio,meanVAF,sdVAF,`5`:`100`) %>% 
        ungroup()
      
      simul_merge <- rbind(simul_merge, temp)
      
    }
    invisible(simul_merge)
  }



simul_result_fixvaf <- setClass("simul_result_fixvaf",
                           slots = c(
                             pseudoVAF = "numeric",
                             VAR = "character",
                             patient = "character",
                             tissue = "character",
                             age = "numeric",
                             fixvaf_simul_original = "tbl_df",
                             fixvaf_simul_original_dist = "tbl_df",
                             fixvaf_simul_sample = "tbl_df",
                             fixvaf_simul_sample_dist = "tbl_df",
                             obs = "tbl_df"
                           ))


file_fixvaf <- getAbsolutePath("~/project/11_Clone_MT/10_Simulation/231220_fe_turnover/","*simul.iter1.tsv")
simul_result_fixvaf_lst <- vector(mode="list", length=length(file_fixvaf))
names(simul_result_fixvaf_lst) <- gsub(".simul.iter1.tsv","",basename(file_fixvaf))


simul_result_fixvaf_list <- mcmapply(function(input){
  
  vafdist_col <- c(seq(5,100,by=5) %>% as.character())
  print(input)
  pat <- unlist(strsplit(basename(input),split="\\_"))[1]
  var_p <- unlist(strsplit(basename(input),split="\\_"))[2]
  tis <- unlist(strsplit(basename(input),split="\\_"))[3]
  ag <- as.numeric(gsub(".simul.iter1.tsv","",unlist(strsplit(basename(input),split="\\_"))[4]))
  pos <- gsub("[ACGTto]","",unlist(strsplit(basename(input),split="\\_"))[2])
  mut <- gsub("to",">",gsub("[0-9]","",unlist(strsplit(basename(input),split="\\_"))[2]))
  var <- paste0(pos," ",mut)
  var_name <- gsub(".simul.iter1.tsv","",basename(input))
  
  # simulation result load
  simul_res_fixvaf <- read_tsv(input) %>% 
    separate(VAFdist,into=vafdist_col, sep=":") %>% 
    mutate(across(`5`:`100`, as.numeric)) %>% 
    mutate(across(`5`:`100`, ~./10000))
  vaf <- simul_res_fixvaf %>% slice(1) %>% pull(initVAF) %>% round(.,5)
  
  # observed_data
  obs <- pharse_obs(df %>% filter(!(sample %in% same_lineage_sample)),pat,var,same_lineage_sample)
  patientN <- sample_info_final_v3 %>% 
    filter(patient==pat) %>% 
    filter(!(sample %in% same_lineage_sample)) %>% 
    filter(type=="clonal_normal") %>% nrow()
  
  # sampling 
  simul_res_sample <- sample_cells_sumstat(df=simul_res_fixvaf, start=5, end=10004, patN = patientN, iter = 100)
  
  # distance
  dist <- calculate_distance(simul_res_fixvaf, obs)
  dist_2 <- calculate_distance(simul_res_sample,obs)
  
  dist <- dist %>% arrange(distance) %>% slice(1:1000)
  dist_2 <- dist_2 %>% arrange(distance) %>% slice(1:1000)
  
  temp_lst <- simul_result_fixvaf(pseudoVAF=vaf,
                           VAR=var,
                           patient=pat,
                           tissue=tis,
                           age=ag,
                           fixvaf_simul_original = simul_res_fixvaf,
                           fixvaf_simul_original_dist = dist,
                           fixvaf_simul_sample = simul_res_sample,
                           fixvaf_simul_sample_dist = dist_2,
                           obs = obs)
  
  simul_result_fixvaf_lst[[var_name]] <- temp_lst
    
  }, input = file_fixvaf, mc.cores=16)

names(simul_result_fixvaf_list) <- gsub(".simul.iter1.tsv","",basename(names(simul_result_fixvaf_list)))
save(simul_result_fixvaf_list, file="~/R/R_script/11_Clone_MT/simul_fixvaf_result_fevar_sample_231221.RData")