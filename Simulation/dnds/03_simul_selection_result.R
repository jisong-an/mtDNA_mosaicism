################################
#
#   Selection simulation result
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

# 1. merge simulation result
# 2. calculate dN/dS of simulated & observed data


# functions
source("~/R/R_script/11_Clone_MT/simul__selection_function.R")
source("~/R/R_script/getAbsolutePath.R")


## 1. analysis -pdf -----------------------------------------------------------


# Analysis::all region

res_file <- getAbsolutePath("~/project/11_Clone_MT/14_Selection/03_simulation/","*allvar.simul_result.tsv")
site_file <- getAbsolutePath("~/project/11_Clone_MT/14_Selection/02_count/", "*tsv")

all_nonsyn_dnds_pval <- tibble(patient=character(), context=character(), dnds_obs = numeric(), upper_pval=numeric(), lower_pval=numeric())
all_nonsense_dnds_pval <- tibble(patient=character(), context=character(), dnds_obs = numeric(), upper_pval=numeric(), lower_pval=numeric())

pdf("~/project/11_Clone_MT/14_Selection/new_230724/all_dnds_simul_result.pdf", width=8, height=5, onefile=T)
for (i in 1:length(res_file)){
  
  res <- read_tsv(res_file_2[i])
  site <- read_tsv(site_file[i])
  pat <- strsplit(tail(strsplit(basename(res_file[i]), split = "_")[[1]],2)[1], "\\.")[[1]][1]
  
  try({
    nonsyn_simul <- parse_simul_result(res, site, "nonsynonymous")
    nonsyn_obs <- make_dnds_observed(pat, "nonsynonymous")
    nonsyn_ctx <- nonsyn_obs %>% pull(context)
    nonsyn_pval <- calculate_pval(nonsyn_simul, nonsyn_obs)
    
    if (nrow(nonsyn_obs)>0){
      
      g1 <- nonsyn_simul %>% 
        filter(context %in% nonsyn_ctx) %>% 
        ggplot() + 
        geom_histogram(aes(x=dnds), bins=15) + 
        geom_vline(nonsyn_obs, mapping=aes(xintercept=dnds), color="darkred")  + 
        facet_wrap(~context, nrow=2) + 
        theme_light() + 
        theme(strip.text.x = element_text(color="black", face="bold", size=10), axis.title.y = element_blank()) + 
        geom_text(nonsyn_pval, mapping=aes(x=10, y=5000,label=paste0("upper_pval : ",pval_upper)), hjust=1, vjust=1) +
        geom_text(nonsyn_pval, mapping=aes(x=10, y=4500,label=paste0("lower_pval : ", pval_lower)), hjust=1, vjust=1) +
        ggtitle(paste0(pat, " : nonsynonymous dN/dS -all region"))
      
      plot(g1)
      
      nonsyn_temp <- nonsyn_obs %>% 
        merge(nonsyn_pval %>% 
                select(context,pval_upper,pval_lower), by=c("context"), all=TRUE) %>% 
        as_tibble() %>% 
        select(context, dnds, pval_upper, pval_lower) %>% 
        rename("dnds_obs"=dnds) %>% 
        mutate(patient=pat)
      
      all_nonsyn_dnds_pval <- rbind(all_nonsyn_dnds_pval, nonsyn_temp)
      
    }
    
  })
  
  try({
    
    nonsense_simul <- parse_simul_result(res,site,"nonsense")
    nonsense_obs <- make_dnds_observed(pat, "nonsense")
    nonsense_ctx <- nonsense_obs %>% pull(context)
    nonsense_pval <- calculate_pval(nonsense_simul, nonsense_obs)
    
    if (nrow(nonsense_obs)>0){
      
      g2 <- nonsense_simul %>% 
        filter(context %in% nonsense_ctx) %>% 
        ggplot() + 
        geom_histogram(aes(x=dnds), bins=15) + 
        geom_vline(nonsense_obs, mapping=aes(xintercept=dnds), color="darkred")  + 
        facet_wrap(~context, nrow=2) + 
        theme_light() + 
        theme(strip.text.x = element_text(color="black", face="bold", size=10), axis.title.y = element_blank()) + 
        geom_text(nonsense_pval, mapping=aes(x=10, y=5000,label=paste0("upper_pval : ",pval_upper)), hjust=1, vjust=1) +
        geom_text(nonsense_pval, mapping=aes(x=10, y=4500,label=paste0("lower_pval : ", pval_lower)), hjust=1, vjust=1) +
        ggtitle(paste0(pat, " : nonsense dN/dS -all region"))
      
      plot(g2)
      
      nonsense_temp <- nonsense_obs %>% 
        merge(nonsense_pval %>% 
                select(context,pval_upper,pval_lower), by=c("context"), all=TRUE) %>% 
        as_tibble() %>% 
        select(context, dnds, pval_upper, pval_lower) %>% 
        rename("dnds_obs"=dnds) %>% 
        mutate(patient=pat)
      
      all_nonsense_dnds_pval <- rbind(all_nonsense_dnds_pval, nonsense_temp)
      
    }
    
  })
  
}
dev.off()



## 2. simul result merge ------------------------------------------------------------

# colon
line1_all_file <- c(getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/03_simulation/","Line1.*allvar.simul_result.tsv"),getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/03_simulation/","Abortus.*allvar.simul_result.tsv"))
site_file <- c(getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/02_count/", "Line1.*tsv"),getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/02_count/", "Abortus.*tsv"))

colon_nonsyn_simul_result_all <- tibble(context=character(), synonymous_n=numeric(), nonsynonymous_n=numeric(), nonsense_n=numeric(), stop_loss_n=numeric(), total=numeric(), synonymous=numeric(), nonsynonymous=numeric(), nonsense=numeric(), stop_loss=numeric(), dnds=numeric())
colon_nonsense_simul_result_all <- tibble(context=character(), synonymous_n=numeric(), nonsynonymous_n=numeric(), nonsense_n=numeric(), stop_loss_n=numeric(), total=numeric(), synonymous=numeric(), nonsynonymous=numeric(), nonsense=numeric(), stop_loss=numeric(), dnds=numeric())


for (i in 1:length(line1_all_file)){
  
  res <- read_tsv(line1_all_file[i])
  site <- read_tsv(site_file[i])
  nonsyn_simul <- parse_simul_result(res, site, "nonsynonymous")
  nonsense_simul <- parse_simul_result(res,site,"nonsense")
  colon_nonsyn_simul_result_all <- rbind(colon_nonsyn_simul_result_all, nonsyn_simul)
  colon_nonsense_simul_result_all <- rbind(colon_nonsense_simul_result_all, nonsense_simul)
  
}


# fibroblast
db_all_file <- getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/03_simulation/","DB.*allvar.simul_result.tsv")
site_file <- c(getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/02_count/", "DB.*tsv"))

fibr_nonsyn_simul_result_all <- tibble(context=character(), synonymous_n=numeric(), nonsynonymous_n=numeric(), nonsense_n=numeric(), stop_loss_n=numeric(), total=numeric(), synonymous=numeric(), nonsynonymous=numeric(), nonsense=numeric(), stop_loss=numeric(), dnds=numeric())
fibr_nonsense_simul_result_all <- tibble(context=character(), synonymous_n=numeric(), nonsynonymous_n=numeric(), nonsense_n=numeric(), stop_loss_n=numeric(), total=numeric(), synonymous=numeric(), nonsynonymous=numeric(), nonsense=numeric(), stop_loss=numeric(), dnds=numeric())


for (i in 1:length(db_all_file)){
  
  res <- read_tsv(db_all_file[i])
  site <- read_tsv(site_file[i])
  nonsyn_simul <- parse_simul_result(res, site, "nonsynonymous")
  nonsense_simul <- parse_simul_result(res,site,"nonsense")
  fibr_nonsyn_simul_result_all <- rbind(fibr_nonsyn_simul_result_all, nonsyn_simul)
  fibr_nonsense_simul_result_all <- rbind(fibr_nonsense_simul_result_all, nonsense_simul)
  
}


# blood
blood_all_file <- c(getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/03_simulation/","Hblood.*allvar.simul_result.tsv"),getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/03_simulation/","blood.*allvar.simul_result.tsv"))
site_file <- c(getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/02_count/", "Hblood.*tsv"),getAbsolutePath("~/project/11_Clone_MT/14_Selection/new_230724/02_count/", "blood.*tsv"))

blood_nonsyn_simul_result_all <- tibble(context=character(), synonymous_n=numeric(), nonsynonymous_n=numeric(), nonsense_n=numeric(), stop_loss_n=numeric(), total=numeric(), synonymous=numeric(), nonsynonymous=numeric(), nonsense=numeric(), stop_loss=numeric(), dnds=numeric())
blood_nonsense_simul_result_all <- tibble(context=character(), synonymous_n=numeric(), nonsynonymous_n=numeric(), nonsense_n=numeric(), stop_loss_n=numeric(), total=numeric(), synonymous=numeric(), nonsynonymous=numeric(), nonsense=numeric(), stop_loss=numeric(), dnds=numeric())

for (i in 1:length(blood_all_file)){
  
  res <- read_tsv(blood_all_file[i])
  site <- read_tsv(site_file[i])
  nonsyn_simul <- parse_simul_result(res, site, "nonsynonymous")
  nonsense_simul <- parse_simul_result(res,site,"nonsense")
  blood_nonsyn_simul_result_all <- rbind(blood_nonsyn_simul_result_all, nonsyn_simul)
  blood_nonsense_simul_result_all <- rbind(blood_nonsense_simul_result_all, nonsense_simul)
  
}



## 3. compare with observation ------------------------------------------

### colon -----------------------------------

nonsyn_obs_colon <- df %>% filter(!is.na(VAR)) %>% 
  filter(project=="Line1" | project=="Abortus") %>% 
  filter(type=="clonal_normal") %>% 
  filter(!(vartype%in%c("fe","fe_gray","lineage"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>%
  select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
  rbind(df %>% filter(project=="Line1" | project=="Abortus") %>% filter(type=="clonal_normal") %>% 
          filter(vartype_2!="frequent_recurrent") %>%
          filter(vartype=="lineage") %>% 
          select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
          unique()) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
  mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
  group_by(context, Consequence) %>% 
  tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n", ifelse(grepl("stop_lost",Consequence),"stop_loss_n",NA))))) %>% 
  select(-Consequence) %>% filter(!is.na(consq))  %>% spread(key="consq", value="n") %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  merge(site %>% rename("context"=VAR), by=c("context")) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(dnds=((nonsynonymous_n/nonsynonymous)/(synonymous_n/synonymous))) %>% 
  filter(!is.infinite(dnds)) %>%
  filter(!is.na(dnds)) %>% 
  filter(dnds!=0) 

nonsense_obs_colon <- df %>% filter(!is.na(VAR)) %>% 
  filter(project=="Line1" | project=="Abortus") %>% 
  filter(type=="clonal_normal") %>% 
  filter(!(vartype%in%c("fe","fe_gray","lineage"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>%
  select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
  rbind(df %>% filter(project=="Line1" | project=="Abortus") %>% filter(type=="clonal_normal") %>% 
          filter(vartype_2!="frequent_recurrent") %>%
          filter(vartype=="lineage") %>% 
          select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
          unique()) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
  mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
  group_by(context, Consequence) %>% 
  tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n", ifelse(grepl("stop_lost",Consequence),"stop_loss_n",NA))))) %>% 
  select(-Consequence) %>% filter(!is.na(consq))  %>% spread(key="consq", value="n") %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  merge(site %>% rename("context"=VAR), by=c("context")) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(dnds=((nonsense_n/nonsense)/(synonymous_n/synonymous))) %>% 
  filter(!is.infinite(dnds)) %>%
  filter(!is.na(dnds)) %>% 
  filter(dnds!=0) 

nonsyn_ctx_colon <- nonsyn_obs_colon %>% filter(nonsynonymous_n > 10) %>% pull(context)
nonsense_ctx_colon <- nonsense_obs_colon %>% filter(nonsense_n > 10) %>% pull(context)



## fibroblast ----------------------------------------------------

nonsyn_obs_fibr <- df %>% filter(!is.na(VAR)) %>% 
  filter(type=="clonal_normal") %>% 
  filter(!(vartype%in%c("fe","fe_gray","lineage"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(project=="DB") %>% 
  select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
  rbind(df %>% filter(project=="DB") %>% filter(!is.na(VAR)) %>% filter(type=="clonal_normal") %>% 
          filter(vartype_2!="frequent_recurrent") %>% 
          filter(vartype=="lineage") %>% 
          select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
          unique()) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
  mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
  group_by(context, Consequence) %>% 
  tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n", ifelse(grepl("stop_lost",Consequence),"stop_loss_n",NA))))) %>% 
  select(-Consequence) %>% filter(!is.na(consq))  %>% spread(key="consq", value="n") %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  merge(site %>% rename("context"=VAR), by=c("context")) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(dnds=((nonsynonymous_n/nonsynonymous)/(synonymous_n/synonymous))) %>% 
  filter(!is.infinite(dnds)) %>%
  filter(!is.na(dnds)) %>% 
  filter(dnds!=0) 

nonsense_obs_fibr <- df %>% filter(!is.na(VAR)) %>% 
  filter(type=="clonal_normal") %>% 
  filter(!(vartype%in%c("fe","fe_gray","lineage"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>% 
  filter(project=="DB") %>% 
  select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
  rbind(df %>% filter(project=="DB") %>% filter(!is.na(VAR)) %>% filter(type=="clonal_normal") %>% 
          filter(vartype_2!="frequent_recurrent") %>% 
          filter(vartype=="lineage") %>% 
          select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
          unique()) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
  mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
  group_by(context, Consequence) %>% 
  tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n", ifelse(grepl("stop_lost",Consequence),"stop_loss_n",NA))))) %>% 
  select(-Consequence) %>% filter(!is.na(consq))  %>% spread(key="consq", value="n") %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  merge(site %>% rename("context"=VAR), by=c("context")) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(dnds=((nonsense_n/nonsense)/(synonymous_n/synonymous))) %>% 
  filter(!is.infinite(dnds)) %>%
  filter(!is.na(dnds)) %>% 
  filter(dnds!=0) 

nonsyn_ctx_fibr <- nonsyn_obs_fibr %>% filter(nonsynonymous_n > 10) %>% filter(context!="A>C") %>% pull(context)
nonsense_ctx_fibr <- nonsense_obs_fibr %>% filter(nonsense_n > 10) %>% pull(context)



## blood ------------------------------------------

nonsyn_obs_blood <- df %>% filter(!is.na(VAR)) %>% 
  filter(project=="Hblood" | project=="blood_new") %>% 
  filter(type=="clonal_normal") %>% 
  filter(!(vartype%in%c("fe","fe_gray","lineage"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>%
  select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
  rbind(df %>% filter(project=="Hblood" | project=="blood_new") %>% filter(type=="clonal_normal") %>% 
          filter(vartype_2!="frequent_recurrent") %>%
          filter(vartype=="lineage") %>% 
          select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
          unique()) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
  mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
  group_by(context, Consequence) %>% 
  tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n", ifelse(grepl("stop_lost",Consequence),"stop_loss_n",NA))))) %>% 
  select(-Consequence) %>% filter(!is.na(consq))  %>% spread(key="consq", value="n") %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  merge(site %>% rename("context"=VAR), by=c("context")) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(dnds=((nonsynonymous_n/nonsynonymous)/(synonymous_n/synonymous))) %>% 
  filter(!is.infinite(dnds)) %>%
  filter(!is.na(dnds)) %>% 
  filter(dnds!=0) 

nonsense_obs_blood <- df %>% filter(!is.na(VAR)) %>% 
  filter(project=="Hblood" | project=="blood_new") %>% 
  filter(type=="clonal_normal") %>% 
  filter(!(vartype%in%c("fe","fe_gray","lineage"))) %>% 
  filter(vartype_2!="frequent_recurrent") %>%
  select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
  rbind(df %>% filter(project=="Hblood" | project=="blood_new") %>% filter(type=="clonal_normal") %>% 
          filter(vartype_2!="frequent_recurrent") %>%
          filter(vartype=="lineage") %>% 
          select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
          unique()) %>% 
  filter(nchar(REF)==nchar(ALT)) %>% 
  mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
  mutate(context=paste0(REF,">",ALT)) %>% filter(consq=="coding")  %>% 
  group_by(context, Consequence) %>% 
  tally() %>% mutate(consq=ifelse(grepl("missense",Consequence), "nonsynonymous_n", ifelse(grepl("synonymous",Consequence), "synonymous_n",ifelse(grepl("stop_gained",Consequence), "nonsense_n", ifelse(grepl("stop_lost",Consequence),"stop_loss_n",NA))))) %>% 
  select(-Consequence) %>% filter(!is.na(consq))  %>% spread(key="consq", value="n") %>% 
  mutate_all(~replace(.,is.na(.),0)) %>% 
  merge(site %>% rename("context"=VAR), by=c("context")) %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(dnds=((nonsense_n/nonsense)/(synonymous_n/synonymous))) %>% 
  filter(!is.infinite(dnds)) %>%
  filter(!is.na(dnds)) %>% 
  filter(dnds!=0) 

nonsyn_ctx_blood <- nonsyn_obs_blood %>% filter(nonsynonymous_n > 15) %>% pull(context)
nonsense_ctx_blood <- nonsense_obs_blood %>% filter(nonsense_n > 10) %>% pull(context)




  
## merge --------------------------

  simul_result_merge_consq <- tibble()
  obs_merge_consq <- tibble()

  for (x in c("colon","fibr","blood")){
    
    nonsyn_simul <- get(paste0(x,"_nonsyn_simul_result_all"))
    nonsense_simul <- get(paste0(x,"_nonsense_simul_result_all"))
    nonsyn_ctx <- get(paste0("nonsyn_ctx_",x))
    nonsense_ctx <- get(paste0("nonsense_ctx_",x))
    nonsyn_obs <- get(paste0("nonsyn_obs_",x))
    nonsense_obs <- get(paste0("nonsense_obs_",x))
    
    temp1 <- nonsyn_simul %>% ungroup() %>% 
      filter(context %in% nonsyn_ctx) %>% 
      arrange(dnds) %>% 
      group_by(context) %>% 
      summarise(q1=quantile(dnds,0.2),q2=quantile(dnds,0.8)) %>% 
      mutate(consequence="nonsynonymous") %>% 
      mutate(tissue=x)
    
    temp2 <- nonsense_simul %>% ungroup() %>% 
      filter(context %in% nonsense_ctx) %>% 
      arrange(dnds) %>% 
      group_by(context) %>% 
      summarise(q1=quantile(dnds,0.2),q2=quantile(dnds,0.8)) %>% 
      mutate(consequence="nonsense") %>% 
      mutate(tissue=x)
    
    temp <- rbind(temp1,temp2)
    simul_result_merge_consq <- rbind(simul_result_merge_consq,temp)
    
    obs <- rbind(nonsyn_obs %>% mutate(consequence="nonsynonymous")%>% mutate(tissue=x), 
                 nonsense_obs %>% mutate(consequence="nonsense")%>% mutate(tissue=x)) 
    obs_merge_consq <- rbind(obs_merge_consq,obs)
    
  }
  
  obs_merge_consq %>% 
    ungroup() %>% 
    filter((nonsynonymous_n > 10 & consequence=="nonsynonymous") | (nonsense_n > 10 & consequence=="nonsense")) %>% 
    filter(context!="A>C") %>% 
    filter(consequence=="nonsynonymous") %>% 
    ggplot()  + 
    geom_errorbar(data=simul_result_merge_consq %>% filter(consequence=="nonsynonymous"), aes(x=context, ymin=q1,ymax=q2), width=0.2, color="gray40") + 
    geom_point(aes(x=context,y=dnds, col=tissue), size=4, shape=18) + 
    theme_classic()+
    facet_wrap(~tissue, scales = "free_x") + 
    scale_color_npg() + 
    ylab("dN/dS") + 
    ggtitle("dN/dS : Nonsynonymous mutations") +
    theme(strip.background = element_rect(fill="lightgray", color=NA), 
          strip.text.x = element_text(face="bold"), title = element_text(size=10))
    
  obs_merge_consq %>% 
    ungroup() %>% 
    filter((nonsynonymous_n > 10 & consequence=="nonsynonymous") | (nonsense_n > 10 & consequence=="nonsense")) %>% 
    filter(context!="A>C") %>% 
    filter(consequence=="nonsense") %>% 
    ggplot()  + 
    geom_errorbar(data=simul_result_merge_consq %>% filter(consequence=="nonsense"), aes(x=context, ymin=q1,ymax=q2), width=0.2, color="gray40") + 
    geom_point(aes(x=context,y=dnds, col=tissue), size=4, shape=18) + 
    theme_classic()+
    facet_wrap(~tissue, scales = "free_x") + 
    scale_color_npg() + 
    ylab("dN/dS") + 
    ggtitle("dN/dS : Nonsense mutations") +
    theme(strip.background = element_rect(fill="lightgray", color=NA), 
          strip.text.x = element_text(face="bold"), title = element_text(size=10))