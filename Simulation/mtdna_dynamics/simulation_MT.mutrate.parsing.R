#######################
#
#  Simulation: mutation rate result parsing
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
#######################

# 1. calculate summary statistics in observed data
# 2. calculate MSE (mean squared error) between simulated & observed data


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
remove_sample <- same_lineage_sample


## simul result : mutation rate  -----------------------------------------------


calculate_distance_vaf <- function(simul,obs){
  
  simul_res <- simul %>%
    mutate(distance=((`5`-obs[["5"]])^2+ 
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
                       (`100`-obs[["100"]])^2+ 
                       (over3N-obs[["over3N"]])^2+ 
                       (var0N-obs[["var0N"]])^2+ 
                       (over2N-obs[["over2N"]])^2)) 
  
  invisible(simul_res)
  
}

simul_result_mrate <- setClass("simul_result_mrate",
                                slots = c(
                                  patient = "character",
                                  tissue = "character",
                                  age = "numeric",
                                  mrate_simul = "data.frame",
                                  mrate_simul_dist = "data.frame",
                                  obs = "data.frame"
                                ))

file_mrate <- getAbsolutePath("~/project/11_Clone_MT/10_Simulation/231226_mrate/merged/","*simul.merged.tsv")
simul_result_mrate_lst <- vector(mode="list", length=length(file_mrate))
names(simul_result_mrate_lst) <- gsub(".simul.merged.tsv","",basename(file_mrate))
vafdist_col <- c(paste0("<",seq(5,100,by=5)),"varnum")
vafrange_df <- as_tibble(matrix(nrow = 0, ncol = length(vafdist_col[-length(vafdist_col)])), 
                         .name_repair = ~ vafdist_col[-length(vafdist_col)])



## simul summary ------------------------------------------------

for (f in file_mrate){
  
  pat <- unlist(strsplit(basename(f),split="\\_"))[1]
  tis <- gsub(".simul.merged.tsv","",unlist(strsplit(basename(f),split="\\_"))[3])
  ag <- as.numeric(unlist(strsplit(basename(f),split="\\_"))[2])
  var_name <- gsub(".simul.merged.tsv","",basename(f))
  
  # simulation result load
  simul_res_mrate_temp <- read_tsv(f, col_names= c("turnover", "age", "mtCN","generation", "sampling","max_cell", "mrate","mrate_gen","result")) %>% 
    mutate(iter=row_number()) %>%                                 # each iteration (1-1000)
    separate_rows(result, sep="__") %>%                           # sample separate
    separate(result, into=c("sample","vafdist"), sep="\\|") %>%   # VAF distribution per sample
    separate(vafdist,into=vafdist_col, sep=":") %>%               # separate each VAF range
    mutate(across((10:ncol(.)), as.numeric)) 
  
  simul_res_mrate_maxvafdist <- simul_res_mrate_temp %>% 
    select(iter,everything()) %>% select(-varnum) %>% 
    gather(11:ncol(.), key="group",value="VARn") %>% 
    group_by(mrate,mrate_gen,iter,sample) %>% 
    arrange(mrate,mrate_gen, iter,sample,factor(group, levels=rev(vafdist_col[-length(vafdist_col)]))) %>% 
    filter(VARn>0) %>% 
    dplyr::slice(1) %>% 
    merge(simul_res_mrate_temp %>% select(turnover:sample,iter), by=c("iter","sample","turnover","age","mtCN","generation","sampling","max_cell","mrate","mrate_gen"), all=TRUE) %>%     # if clone has no variant, assign to "<5"
    as_tibble() %>%
    mutate(across(group, ~replace(.,is.na(.),"<5"))) %>% 
    mutate(across(VARn, ~replace(.,is.na(.),1))) %>% 
    group_by(mrate,mrate_gen,iter,group) %>% 
    tally() %>% 
    arrange(mrate,mrate_gen,iter,factor(group, levels=vafdist_col[-length(vafdist_col)])) %>% 
    spread(key="group",value="n") %>% 
    merge(vafrange_df, all=TRUE) %>% 
    as_tibble() %>% 
    arrange(iter) %>% 
    select(mrate,mrate_gen, vafdist_col[-length(vafdist_col)]) %>% 
    mutate_all(~replace(.,is.na(.),0))
  
  simul_res_mrate_homoprop <- simul_res_mrate_temp %>% 
    select(mrate,mrate_gen,iter,sample,`<100`,varnum) %>% 
    group_by(mrate,mrate_gen, iter) %>% 
    mutate(over2=rowSums(across(`<100`)>=2)) %>% 
    mutate(over3=rowSums(across(`<100`)>=3)) %>% 
    mutate(var0=rowSums(across(varnum)==0)) %>% 
    mutate(over2N=sum(over2), over3N=sum(over3), var0N=sum(var0)) %>% 
    select(mrate,mrate_gen,iter,over2N,over3N,var0N) %>% 
    unique() %>% 
    ungroup() %>% 
    select(over2N,over3N,var0N)
  
  simul_res_mrate <- cbind(simul_res_mrate_maxvafdist, simul_res_mrate_homoprop) %>% as_tibble() %>%
    rename_at(3:22,~gsub("<","",.))
  
  
  # observed_data
  maxvaf_dist_obs <- df %>% 
    filter(patient == pat) %>% 
    filter(!(sample %in% remove_sample)) %>% 
    filter(type=="clonal_normal") %>% 
    filter(vartype!="fe") %>% 
    filter(vartype_2 !=" frequent_recurrent") %>% 
    filter(POS!=414) %>% 
    filter(!is.na(VAF)) %>% 
    group_by(patient,sample) %>% 
    summarise(max=max(VAF,na.rm=TRUE)) %>% 
    merge(sample_info_final_v3 %>% 
            filter(patient == pat) %>% 
            filter(!(sample %in% remove_sample)) %>% 
            filter(type=="clonal_normal") %>% 
            select(patient,sample), by=c("patient","sample"), all=TRUE) %>% 
    as_tibble() %>% 
    mutate(across(max, ~replace(.,is.na(.),1))) %>% 
    mutate(x=(max%/%5)+1) %>% 
    mutate(group=paste0("<",x*5)) %>% 
    group_by(patient,group) %>% 
    tally() %>% 
    arrange(patient,factor(group, levels=vafdist_col[-length(vafdist_col)])) %>% 
    spread(key="group",value="n") %>% 
    merge(vafrange_df, all=TRUE) %>% 
    as_tibble() %>% 
    select(patient, vafdist_col[-length(vafdist_col)]) %>% 
    mutate_all(~replace(.,is.na(.),0)) %>% 
    dplyr::rename("mrate_gen"=patient) %>% arrange(mrate_gen)
  
  homoprop_obs <- df %>% 
    filter(patient == pat) %>% 
    filter(!(sample %in% remove_sample)) %>% 
    filter(type=="clonal_normal") %>% 
    filter(vartype!="fe") %>% 
    filter(vartype_2 !=" frequent_recurrent") %>% 
    filter(POS!=414) %>% 
    filter(!is.na(VAF)) %>% 
    group_by(patient,sample) %>% 
    filter(VAF>90) %>% 
    tally() %>% 
    group_by(patient) %>% 
    mutate(over2=rowSums(across(n)>=2)) %>% 
    mutate(over3=rowSums(across(n)>=3)) %>% 
    mutate(over2N=sum(over2),over3N=sum(over3)) %>% 
    select(patient,over2N,over3N) %>% unique() %>% ungroup() %>% 
    merge(tibble(patient=pat),all=TRUE) %>% 
    as_tibble() %>% 
    mutate(across(over2N:over3N, ~replace(.,is.na(.),0))) %>% select(over2N,over3N)
  
  var0_obs <- df %>% 
    filter(patient == pat) %>% 
    filter(!(sample %in% remove_sample)) %>% 
    filter(type=="clonal_normal") %>% 
    filter(vartype!="fe") %>% 
    filter(vartype_2 !=" frequent_recurrent") %>% 
    filter(POS!=414)  %>% 
    merge(sample_info_final_v3 %>% 
            filter(patient == pat) %>% 
            filter(type=="clonal_normal") %>% 
            select(patient,sample), by=c("patient","sample"), all=TRUE) %>% 
    as_tibble() %>% 
    filter(is.na(VAR)) %>% 
    group_by(patient) %>% 
    tally() %>% 
    dplyr::rename("var0N"=n) %>% 
    merge(tibble(patient=pat), all=TRUE) %>% 
    as_tibble() %>% 
    mutate(across(var0N, ~replace(.,is.na(.),0))) %>% select(var0N)
  

  res_obs <- cbind(maxvaf_dist_obs,homoprop_obs) %>% cbind(var0_obs)%>% as_tibble() %>% 
    rename_at(2:21,~gsub("<","",.))

  # distance
  dist <- calculate_distance_vaf(simul_res_mrate, res_obs)
  dist <- dist %>% arrange(distance) %>% slice(1:1000)
  
  temp_lst <- simul_result_mrate( patient=pat,
                                  tissue=tis,
                                  age=ag,
                                  mrate_simul = simul_res_mrate,
                                  mrate_simul_dist = dist,
                                  obs = res_obs)
  
  simul_result_mrate_lst[[var_name]] <- temp_lst

}




mrate_dist_merge <- do.call(rbind, mclapply(simul_result_mrate_lst, function(x) {
  x@mrate_simul_dist <- x@mrate_simul_dist %>% dplyr::slice(1:50) %>% 
    mutate(patient=x@patient, age=x@age, tissue=x@tissue)
  return(x@mrate_simul_dist)
}, mc.cores=2))


save(simul_result_mrate_lst, mrate_dist_merge, file="~/R/R_script/11_Clone_MT/simul_result_mrate_231229.RData")