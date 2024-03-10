################################
#
#   make context list
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

## make observed context count list per individual

pat_list <- sample_info_final_v3 %>% filter(project%in%c("Line1","DB","Hblood", "blood_new","Abortus")) %>% select(patient) %>% unique() %>% pull(patient)
context_list=c("A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G")

for (pat in pat_list){
  
  proj <- sample_info_final_v3 %>% filter(patient==pat) %>% pull(project) %>% .[1]
  
  df %>% filter(type %in% c("clonal_normal")) %>% 
    filter(patient==pat) %>% 
    filter(!(vartype%in%c("fe","fe_gray","lineage")))%>% 
    filter(!is.na(VAR)) %>% 
    filter(vartype_2!="frequent_recurrent") %>% 
    select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
    rbind(df %>% filter(patient==pat) %>% 
            filter(vartype=="lineage") %>% 
            select(project,patient,VAR:Gene,gender:tissue, vartype) %>% 
            unique()) %>% 
    filter(nchar(REF)==nchar(ALT)) %>% 
    mutate(consq=ifelse(grepl("gene",Consequence), "noncoding","coding")) %>% 
    mutate(context=paste0(REF,">",ALT)) %>% 
    group_by(context) %>% 
    tally() %>% 
    merge(tibble(context=context_list), by=c("context"), all=TRUE) %>% 
    as_tibble() %>% 
    mutate_all(~replace(., is.na(.),0)) %>% 
    write_tsv(paste0("~/project/11_Clone_MT/14_Selection/new_230724/03_simulation/", proj,"_",pat,".allvar.tsv"))
  
}
