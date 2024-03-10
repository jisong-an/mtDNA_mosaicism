################################
#
#   Functions for analysis of selection simulation
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################


# calculate dnds from observed data
make_dnds_observed <- function(pat, mut_consq){
  
  df_temp <- df %>% filter(patient==pat) %>% filter(type %in% c("clonal_normal","adenoma")) %>% 
    filter(!(vartype%in%c("fe", "fe_gray","lineage"))) %>% 
    filter(!is.na(VAR)) %>% 
    filter(vartype_2!="frequent_recurrent") %>% 
    select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
    rbind(df %>% filter(patient==pat) %>% 
            filter(vartype=="lineage") %>% 
            select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
            unique()) %>% 
    filter(nchar(REF)==nchar(ALT)) %>% 
    mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
    mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
    group_by(context, Consequence) %>% 
    tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n","stop_loss_n")))) %>% 
    select(-Consequence) %>% 
    spread(key="consq", value="n") %>% 
    mutate_all(~replace(.,is.na(.),0)) %>% 
    merge(site %>% rename("context"=VAR), by=c("context")) %>% 
    as_tibble() %>% 
    rowwise() %>% 
    mutate(dnds=((get(paste0(mut_consq,"_n"))/get(mut_consq))/(synonymous_n/synonymous))) %>% 
    filter(!is.infinite(dnds)) %>%
    filter(!is.na(dnds)) %>% 
    filter(dnds!=0) 
  
  invisible(df_temp)
  
}


# parse simulation result & calculate dnds
parse_simul_result <- function(df, site, mut_consq){
  
  context_list=c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")
  
  df_parsed <- df %>% 
    separate_rows(synonymous_var, sep="\\|") %>% 
    separate(synonymous_var, sep=":", into=c("context","synonymous_n")) %>% 
    select(iter:synonymous_n) %>% 
    merge(df %>% 
            separate_rows(nonsynonymous_var, sep="\\|") %>% 
            separate(nonsynonymous_var, sep=":", into=c("context","nonsynonymous_n")) %>% 
            select(iter:stop_loss,context,nonsynonymous_n), by=c("iter","total","synonymous","nonsynonymous","nonsense","stop_loss","context"), all=TRUE) %>% 
    as_tibble() %>% 
    merge(df %>% 
            separate_rows(nonsense_var, sep="\\|") %>% 
            separate(nonsense_var, sep=":", into=c("context","nonsense_n")) %>% 
            select(iter:stop_loss,context,nonsense_n), by=c("iter","total","synonymous","nonsynonymous","nonsense","stop_loss","context"), all=TRUE) %>% 
    as_tibble() %>% 
    merge(df %>% 
            separate_rows(stop_loss_var, sep="\\|") %>% 
            separate(stop_loss_var, sep=":", into=c("context","stop_loss_n")) %>% 
            select(iter:stop_loss,context,stop_loss_n), by=c("iter","total","synonymous","nonsynonymous","nonsense","stop_loss","context"), all=TRUE) %>% 
    as_tibble() %>% 
    filter(!is.na(context)) %>% 
    merge(tibble(context=context_list), by=c("context"), all=TRUE) %>% 
    as_tibble() %>% 
    mutate(across(c(synonymous_n:stop_loss_n), as.numeric)) %>% 
    mutate_all(~replace(., is.na(.),0)) %>% 
    select(context, synonymous_n:stop_loss_n) %>% 
    merge(site %>% rename("context"=VAR), by=c("context")) %>% 
    as_tibble() %>% 
    rowwise() %>% 
    mutate(dnds=((get(paste0(mut_consq,"_n"))/get(mut_consq))/(synonymous_n/synonymous))) %>% 
    filter(!is.infinite(dnds)) %>% 
    filter(!is.na(dnds))
  
  invisible(df_parsed)
  
}



# calculate p-value
calculate_pval <- function(simul_res, obs_res){
  
  pval_upper <- simul_res %>% 
    merge(obs_res %>% 
            select(context,dnds) %>% 
            rename("dnds_obs"=dnds), by=c("context")) %>% 
    as_tibble() %>% 
    filter(dnds>=dnds_obs) %>% 
    group_by(context) %>% 
    tally() %>% 
    mutate(pval_upper=n/10000)
  
  pval_lower <- simul_res %>% 
    merge(obs_res %>% 
            select(context,dnds) %>% 
            rename("dnds_obs"=dnds), by=c("context")) %>% 
    as_tibble() %>% 
    filter(dnds<=dnds_obs) %>% 
    group_by(context) %>% 
    tally() %>% 
    mutate(pval_lower=n/10000)
  
  pval <- pval_upper %>% merge(pval_lower, by=c("context"), all=TRUE) %>% as_tibble()
  
  return(pval)
  
}



