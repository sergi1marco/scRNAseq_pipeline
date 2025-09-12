library(tidyverse)
library(rstatix)
library(ggalluvial)
library(gt)
knitr::opts_chunk$set(dev = c("png","pdf"), fig.path = "./Figs/", fig.retina = 3)
rawDf <- read_tsv("Data/Raw/230920_batch1_multiplex_data.tsv") 

tidiedDf <- rawDf %>% 
  select(ImageName = Name, 
         FOV= Image,
         Phenotype = `ObjectInfo - LabelName`,
         Region = `ObjectInfo - ROIName`,
         x = CellX,
         y= CellY,
  ) %>% 
  filter(!str_detect(FOV,"027.+PLEX_5|055.+PLEX_1|046.+53140")) %>% 
  filter(Phenotype != "Class 2",
         Region != "Clear") %>% 
  mutate(ImageName = str_replace(ImageName,
                                 "LIVER_RB045-UB000651",
                                 "LIVER_RB045-2-UB000651"),
         ImageName= str_replace(ImageName, "(RB[0-9]{2,3})(-[A-Z])(.+)", "\\1-1\\2"),
         PatientID = str_extract(ImageName, "RB[0-9]{2,3}-[12]"),
         PatientID = str_replace(PatientID, "(RB)([0-9]{2})-", "\\10\\2-"),
         Phenotype = str_remove_all(Phenotype, "\\(Opal [0-9]{3}\\)| "),
         FOV = str_replace_all(FOV,"\\\\","/") %>%  basename())          
total<-tidiedDf %>% 
  group_by(PatientID,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  mutate(Region = "total") %>% 
  ungroup()

regional<-tidiedDf %>% 
  group_by(PatientID,Region,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  ungroup()
combined <- rbind(regional,total) %>% arrange(PatientID)

combined %>% write_csv("Data/combinedSummaryBatch1Data.csv")
combined %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma")) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1))+
  theme_minimal()+
  facet_wrap(vars(name), scales = "free_y")
rawDf <- read_tsv("Data/Raw/230921_batch2_multiplex_data.tsv", col_select = c(ImageName = Name, 
                                                                              FOV= Image,
                                                                              Phenotype = `ObjectInfo - LabelName`,
                                                                              Region = `ObjectInfo - ROIName`,
                                                                              x = CellX,
                                                                              y= CellY,
)) 

tidiedDf <- rawDf %>% 
  filter(!str_detect(FOV,"1042.+23038|660931.+20912")) %>% 
  filter(Phenotype != "Class 2",
         Region != "Clear") %>% 
  mutate(
    PatientID = str_extract(ImageName,"LIVER_([0-9]+)-", group =1),
    Phenotype = str_remove_all(Phenotype, "\\(Opal [0-9]{3}\\)| ")
  ) 

total<-tidiedDf %>% 
  group_by(PatientID,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  mutate(Region = "total") %>% 
  ungroup()

regional<-tidiedDf %>% 
  group_by(PatientID,Region,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  ungroup()
combined <- rbind(regional,total) %>% arrange(PatientID)
combined %>% write_csv("Data/combinedSummaryBatch2Data.csv")

allSamples <- rbind(read_csv("Data/combinedSummaryBatch1Data.csv"),
                    read_csv("Data/combinedSummaryBatch2Data.csv")) %>% 
  mutate(PatientID = case_when(PatientID == "61434" ~ "614634",
                               PatientID == "59640" ~ "596640",
                               PatientID == "660931" ~ "66093",
                               .default = PatientID))

cohortData <- read_csv("Data/CohortData/subtype_and_sample_type.csv") %>% 
  mutate(PatientID = ifelse(str_detect(PatientID,"RB") & !str_detect(PatientID,"-"),
                            paste0(PatientID,"-1"),
                            PatientID))

allSamples<-allSamples %>%  left_join(cohortData) %>% write_csv("Data/allSamples.csv")
allSamples<-read_csv("Data/allSamples.csv")

allSamples %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()

allSamples %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name) %>% 
  wilcox_test(value~Region, paired = T) %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif  ~ "") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples<-read_csv("Data/allSamples.csv")
allSamples %>%
  select(PatientID, Region,contains("percent"), `sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA, aes(color = `sample type`))+
  geom_point(position = position_jitterdodge(), shape =16, size = 0.5, aes(color = `sample type`))+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()

allSamples %>%
  select(PatientID, Region,contains("percent"), st =`sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name, Region) %>% 
  wilcox_test(value~st, paired = F) %>%
  add_significance() %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             Region ~ "",
             p.signif ~ "") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1, Region))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma")) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  facet_wrap(vars(name), scales = "free_y")

allSamples %>%
  select(PatientID,subtype, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma")) %>%  
  ggplot(aes(subtype,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  facet_wrap(vars(name, Region), scales = "free_x")+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  coord_flip()

allSamples %>%
  select(PatientID,subtype, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region=="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         subtype = str_replace(subtype," ", "\n"),
         subtype = fct_relevel(subtype,"iCCA", "pCCA","dCCA","mixed", "GBC"),
         subtype = fct_rev(subtype)) %>%  
  ggplot(aes(subtype,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  facet_wrap(vars(name), scales = "free_x")+
  coord_flip()+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  labs(y= NULL,x = NULL)

allSamples %>%
  select(PatientID, Region,contains("percent"), subtype) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name) %>% 
  wilcox_test(value~subtype, paired = F) %>%
  add_significance() %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif~"") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples %>% 
  nest(.by = Region) %>% 
  pwalk(function(Region,data){
    
    p<-data %>%
      select(PatientID,subtype,contains("percent")) %>%
      pivot_longer(cols = contains("percent")) %>% 
      filter(!str_detect(name,"Uni|n_")) %>% 
      mutate(name = str_remove(name,"_percent"),
             subtype = str_replace(subtype," ", "\n"),
             subtype = fct_relevel(subtype,"iCCA", "pCCA","dCCA","mixed", "GBC"),
             subtype = fct_rev(subtype)) %>%  
      ggplot(aes(subtype,value))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
      theme_minimal()+
      facet_wrap(vars(name), scales = "free_x")+
      coord_flip()+
      theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"),
            plot.title = element_text(face = "bold"),
            plot.title.position = "plot")+
      labs(title = str_to_title(Region),y= NULL,x = NULL)
    print(p)
  }
  )

allSamples %>%
  select(PatientID, Region,contains("percent"), subtype) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  filter(!str_detect(name,"FAP") ) %>%
  group_by(name, Region, subtype) %>% 
  mutate(sumZero = sum(value)) %>% 
  filter(sumZero!=0) %>% 
  group_by(name, Region) %>%
  wilcox_test(value~subtype, paired = F) %>%
  add_significance() %>% 
  group_by(Region, name) %>% 
  gt() %>% 
  cols_hide(c(.y.,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif~"") %>% 
  # fmt_scientific(columns = -c(n1,n2)) %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples<-read_csv("Data/allSamples.csv") %>% filter(str_detect(subtype, "CCA"))
allSamples %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()

allSamples %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name) %>% 
  wilcox_test(value~Region, paired = T) %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif  ~ "") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples %>%
  select(PatientID, Region,contains("percent"), `sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA, aes(color = `sample type`))+
  geom_point(position = position_jitterdodge(), shape =16, size = 0.5, aes(color = `sample type`))+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()

allSamples %>%
  select(PatientID, Region,contains("percent"), st =`sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name, Region) %>% 
  wilcox_test(value~st, paired = F) %>%
  add_significance() %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             Region ~ "",
             p.signif ~ "") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1, Region))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples %>%
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma")) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  facet_wrap(vars(name), scales = "free_y")

allSamples %>%
  select(PatientID,subtype, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma")) %>%  
  ggplot(aes(subtype,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  facet_wrap(vars(name, Region), scales = "free_x")+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  coord_flip()

allSamples %>%
  select(PatientID,subtype, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region=="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         subtype = str_replace(subtype," ", "\n"),
         subtype = fct_relevel(subtype,"iCCA", "pCCA","dCCA","mixed", "GBC"),
         subtype = fct_rev(subtype)) %>%  
  ggplot(aes(subtype,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
  theme_minimal()+
  facet_wrap(vars(name), scales = "free_x")+
  coord_flip()+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  labs(y= NULL,x = NULL)

allSamples %>%
  select(PatientID, Region,contains("percent"), subtype) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name) %>% 
  wilcox_test(value~subtype, paired = F) %>%
  add_significance() %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif~"") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

allSamples %>% 
  filter(str_detect(subtype,"CCA")) %>% 
  nest(.by = Region) %>% 
  pwalk(function(Region,data){
    
    p<-data %>%
      select(PatientID,subtype,contains("percent")) %>%
      pivot_longer(cols = contains("percent")) %>% 
      filter(!str_detect(name,"Uni|n_")) %>% 
      mutate(name = str_remove(name,"_percent"),
             subtype = str_replace(subtype," ", "\n"),
             subtype = fct_relevel(subtype,"iCCA", "pCCA","dCCA","mixed", "GBC"),
             subtype = fct_rev(subtype)) %>%  
      ggplot(aes(subtype,value))+
      geom_boxplot(outlier.shape = NA)+
      geom_point(position = position_jitter(width = 0.1), shape =16, size = 0.5)+
      theme_minimal()+
      facet_wrap(vars(name), scales = "free_x")+
      coord_flip()+
      theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"),
            plot.title = element_text(face = "bold"),
            plot.title.position = "plot")+
      labs(title = str_to_title(Region),y= NULL,x = NULL)
    print(p)
  }
  )

allSamples %>%
  select(PatientID, Region,contains("percent"), subtype) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         str_detect(subtype,"CCA"),
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  filter(!str_detect(name,"FAP") ) %>%
  group_by(name, Region, subtype) %>% 
  mutate(sumZero = sum(value)) %>% 
  filter(sumZero!=0) %>% 
  group_by(name, Region) %>%
  wilcox_test(value~subtype, paired = F) %>%
  add_significance() %>% 
  group_by(Region, name) %>% 
  gt() %>% 
  cols_hide(c(.y.,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif~"") %>% 
  # fmt_scientific(columns = -c(n1,n2)) %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

library(dplyr)
library(ggplot2)

df <- allSamples %>% 
  select(PatientID, Region,contains("percent"),
         -matches("n_|Uni")) %>% 
  filter(!str_detect(Region,"total")) %>% 
  rename_all(~str_remove(.,"_percent")) %>% 
  #pivot_wider(names_from = Region, id_cols = PatientID, values_from = where(is.numeric)) %>% 
  pivot_longer(cols = where(is.numeric), names_to = c("Phenotype")) %>%
  group_by(Region, Phenotype) %>% 
  summarise(across(where(is.numeric),mean)) %>% 
  ungroup() %>% 
  mutate(Region = str_to_title(Region) %>% str_replace(" ","\n"),
         Region = fct_relevel(Region,"Tumour", "Intratumoural\nStroma","Peritumoural\nStroma","Distal\nStroma"),
         Phenotype = fct_relevel(Phenotype))

orderedLevels <- df %>% filter(str_detect(Region,"Tum")) %>% mutate(Phenotype = fct_reorder(Phenotype, value) %>% fct_rev()) %>% pull(Phenotype) %>% levels()

df %>% 
  mutate(Phenotype = fct_relevel(Phenotype, orderedLevels)) %>% 
  ggplot(aes(x= Region, y = value, alluvium = Phenotype))+
  geom_alluvium(aes(fill = Phenotype), alpha =0.75, decreasing = F, width = 0.01)+
  theme_void()+
  theme(axis.text.x = element_text(margin = margin(t =-20)),
        plot.margin = margin(b = 10))

df <- allSamples %>% 
  select(PatientID, Region,contains("percent"),
         -matches("n_|Uni")) %>% 
  filter(!str_detect(Region,"total")) %>% 
  rename_all(~str_remove(.,"_percent")) %>% 
  pivot_longer(cols = where(is.numeric), names_to = c("Phenotype")) %>%
  group_by(Region, Phenotype) %>% 
  summarise(across(where(is.numeric),mean)) %>%
  filter(str_detect(Phenotype,"CD")) %>% 
  ungroup() %>% 
  mutate(Region = str_to_title(Region) %>% str_replace(" ","\n"),
         Region = fct_relevel(Region,"Tumour", "Intratumoural\nStroma","Peritumoural\nStroma","Distal\nStroma"),
         Phenotype = fct_relevel(Phenotype),
         Phenotype = fct_reorder(Phenotype, value) %>% fct_rev())

orderedLevels <- df %>% filter(str_detect(Region,"Dist")) %>% mutate(Phenotype = fct_reorder(Phenotype, value) %>% fct_rev()) %>% pull(Phenotype) %>% levels()  

df %>% 
  mutate(Phenotype = fct_relevel(Phenotype, orderedLevels)) %>% 
  ggplot(aes(x= Region, y = value, alluvium = Phenotype))+
  geom_alluvium(aes(fill = Phenotype), alpha =0.75, decreasing = F, width = 0.02)+
  theme_void()+
  theme(axis.text.x = element_text(),
        plot.margin = margin(b = 10))

df <- allSamples %>% 
  select(PatientID, Region,contains("percent"),
         -matches("n_|Uni")) %>% 
  filter(!str_detect(Region,"total")) %>% 
  rename_all(~str_remove(.,"_percent")) %>% 
  pivot_longer(cols = where(is.numeric), names_to = c("Phenotype")) %>%
  group_by(Region, Phenotype) %>% 
  summarise(across(where(is.numeric),mean)) %>%
  filter(str_detect(Phenotype,"FAP|SMA|^FSP$")) %>% 
  ungroup() %>% 
  mutate(Region = str_to_title(Region) %>% str_replace(" ","\n"),
         Region = fct_relevel(Region,"Tumour", "Intratumoural\nStroma","Peritumoural\nStroma","Distal\nStroma"),
         Phenotype = fct_relevel(Phenotype),
         Phenotype = fct_reorder(Phenotype, value) %>% fct_rev())

orderedLevels <- df %>% filter(str_detect(Region,"Dist")) %>% mutate(Phenotype = fct_reorder(Phenotype, value) %>% fct_rev()) %>% pull(Phenotype) %>% levels()  

df %>% 
  mutate(Phenotype = fct_relevel(Phenotype, orderedLevels)) %>% 
  ggplot(aes(x= Region, y = value, alluvium = Phenotype))+
  geom_alluvium(aes(fill = Phenotype), alpha =0.75, decreasing = F, width = 0.02)+
  theme_void()+
  theme(axis.text.x = element_text(),
        plot.margin = margin(b = 10))

Batch1Data <- read_tsv("Data/Raw/230920_batch1_multiplex_data.tsv")%>% 
  select(ImageName = Name, 
         FOV= Image,
         Phenotype = `ObjectInfo - LabelName`,
         Region = `ObjectInfo - ROIName`,
         x = CellX,
         y= CellY,
  ) %>% 
  filter(!str_detect(FOV,"027.+PLEX_5|055.+PLEX_1|046.+53140")) %>% 
  filter(Phenotype != "Class 2",
         Region != "Clear") %>% 
  mutate(ImageName = str_replace(ImageName,
                                 "LIVER_RB045-UB000651",
                                 "LIVER_RB045-2-UB000651"),
         ImageName= str_replace(ImageName, "(RB[0-9]{2,3})(-[A-Z])(.+)", "\\1-1\\2"),
         PatientID = str_extract(ImageName, "RB[0-9]{2,3}-[12]"),
         PatientID = str_replace(PatientID, "(RB)([0-9]{2})-", "\\10\\2-"),
         Phenotype = str_remove_all(Phenotype, "\\(Opal [0-9]{3}\\)| "),
         FOV = str_replace_all(FOV,"\\\\","/") %>%  basename())          

Batch2Data <- read_tsv("Data/Raw/230921_batch2_multiplex_data.tsv", col_select = c(ImageName = Name, 
                                                                                   FOV= Image,
                                                                                   Phenotype = `ObjectInfo - LabelName`,
                                                                                   Region = `ObjectInfo - ROIName`,
                                                                                   x = CellX,
                                                                                   y= CellY,
))%>% 
  filter(!str_detect(FOV,"1042.+23038|660931.+20912")) %>% 
  filter(Phenotype != "Class 2",
         Region != "Clear") %>% 
  mutate(
    PatientID = str_extract(ImageName,"LIVER_([0-9]+)-", group =1),
    Phenotype = str_remove_all(Phenotype, "\\(Opal [0-9]{3}\\)| ")
  ) 

allCells <- rbind(Batch1Data,Batch2Data) %>% 
  mutate(PatientID = case_when(PatientID == "61434" ~ "614634",
                               PatientID == "59640" ~ "596640",
                               PatientID == "660931" ~ "66093",
                               .default = PatientID))
rm("Batch1Data","Batch2Data")
gc()

allCells<- allCells %>% 
  mutate(RefinedRegion = ifelse((str_detect(Region, "^tumour") &
                                   Phenotype != "CK"),
                                "intratumoural stroma",Region))

cohortData <- read_csv("Data/CohortData/subtype_and_sample_type.csv") %>% 
  mutate(PatientID = ifelse(str_detect(PatientID,"RB") & !str_detect(PatientID,"-"),
                            paste0(PatientID,"-1"),
                            PatientID))

allCells<-allCells %>%  left_join(cohortData) 
allCells %>% write_rds("Data/allCells.rds")

allCells <- read_rds("Data/allCells.rds")
allCells <- allCells %>% 
  mutate(RefinedRegion = ifelse((str_detect(Region, "^tumour") &
                                   Phenotype != "CK"),
                                "intratumoural stroma",Region),
         RefinedRegion = ifelse((str_detect(Region,"intra") &
                                   Phenotype == "CK"),
                                "tumour", RefinedRegion))

total<-allCells %>% 
  group_by(PatientID,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  mutate(RefinedRegion = "total") %>% 
  ungroup()

regional<-allCells %>% 
  group_by(PatientID,RefinedRegion,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  ungroup()
combined <- rbind(regional,total) %>% arrange(PatientID)

cohortData <- read_csv("Data/CohortData/subtype_and_sample_type.csv") %>% 
  mutate(PatientID = ifelse(str_detect(PatientID,"RB") & !str_detect(PatientID,"-"),
                            paste0(PatientID,"-1"),
                            PatientID))

combined<-combined %>%  left_join(cohortData)

combined %>%
  select(PatientID, RefinedRegion,contains("percent"), `sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(RefinedRegion!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         RefinedRegion = str_replace(RefinedRegion," ", "\n"),
         RefinedRegion = fct_relevel(RefinedRegion,"tumour", "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         RefinedRegion = fct_rev(RefinedRegion),
         RefinedRegion=fct_relabel(RefinedRegion,str_to_title)) %>%  
  ggplot(aes(RefinedRegion,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), shape =16, size = 0.5)+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()

allCells <- read_rds("Data/allCells.rds") %>% select(-RefinedRegion)

total<-allCells %>% 
  filter(Phenotype != "CK", Region != "tumour") %>% 
  group_by(PatientID,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  mutate(Region = "total") %>% 
  ungroup()

regional<-allCells %>% 
  filter(Phenotype != "CK", Region != "tumour") %>% 
  group_by(PatientID,Region,Phenotype) %>% 
  count() %>% 
  pivot_wider(values_from = n, names_from = Phenotype, values_fill = 0) %>% 
  mutate(n= sum(across(where(is.numeric)))) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric),
                ~100*.x/n,
                .names = "{.col}_percent"
  )) %>% 
  ungroup()
combined <- rbind(regional,total) %>% arrange(PatientID)

cohortData <- read_csv("Data/CohortData/subtype_and_sample_type.csv") %>% 
  mutate(PatientID = ifelse(str_detect(PatientID,"RB") & !str_detect(PatientID,"-"),
                            paste0(PatientID,"-1"),
                            PatientID))

combined<-combined %>%  left_join(cohortData)

combined %>%
  select(PatientID, Region,contains("percent"), `sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region, "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(), shape =16, size = 0.5)+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()
combined %>% 
  select(PatientID, Region,contains("percent")) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region, "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name) %>% 
  wilcox_test(value~Region, paired = T) %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             p.adj.signif  ~ "") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.adj.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

combined %>%
  select(PatientID, Region,contains("percent"), `sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region, "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>%  
  ggplot(aes(Region,value))+
  geom_boxplot(outlier.shape = NA, aes(color = `sample type`))+
  geom_point(position = position_jitterdodge(), shape =16, size = 0.5, aes(color = `sample type`))+
  theme_minimal()+
  labs(y=NULL,x=NULL)+
  theme(strip.text = element_text(face = "bold"),panel.spacing.y = unit(10,"mm"))+
  facet_wrap(vars(name), scales = "free_x")+coord_flip()

combined %>%
  select(PatientID, Region,contains("percent"), st =`sample type`) %>%
  pivot_longer(cols = contains("percent")) %>% 
  filter(Region!="total",
         !str_detect(name,"Uni|n_")) %>% 
  mutate(name = str_remove(name,"_percent"),
         Region = str_replace(Region," ", "\n"),
         Region = fct_relevel(Region, "intratumoural\nstroma","peritumoural\nstroma","distal\nstroma"),
         Region = fct_rev(Region),
         Region=fct_relabel(Region,str_to_title)) %>% 
  group_by(name, Region) %>% 
  wilcox_test(value~st, paired = F) %>%
  add_significance() %>% 
  group_by(name) %>% 
  gt() %>% 
  cols_hide(c(.y.,n1,n2,statistic)) %>%
  cols_merge(c(group1,group2), pattern = "{1} vs {2}") %>% 
  cols_label(group1 ~ "",
             Region ~ "",
             p.signif ~ "") %>% 
  fmt_scientific() %>% 
  tab_style(cell_text(weight = "bold"), cells_group()) %>%
  tab_style(cell_text(align = "left"), cells_body(columns = c(group1, Region))) %>% 
  tab_style(cell_fill(color = "#9BCD9B", alpha = 0.5), cells_body(rows = !str_detect(p.signif,"ns"))) %>% 
  tab_options(table.border.top.style = "hidden") 

df <- combined %>% 
  select(PatientID, Region,contains("percent"),
         -matches("n_|Uni")) %>% 
  filter(!str_detect(Region,"total")) %>% 
  rename_all(~str_remove(.,"_percent")) %>% 
  pivot_longer(cols = where(is.numeric), names_to = c("Phenotype")) %>%
  group_by(Region, Phenotype) %>% 
  summarise(across(where(is.numeric),mean)) %>%
  filter(str_detect(Phenotype,"CD")) %>% 
  ungroup() %>% 
  mutate(Region = str_to_title(Region) %>% str_replace(" ","\n"),
         Region = fct_relevel(Region,"Tumour", "Intratumoural\nStroma","Peritumoural\nStroma","Distal\nStroma"),
         Phenotype = fct_relevel(Phenotype),
         Phenotype = fct_reorder(Phenotype, value) %>% fct_rev())

orderedLevels <- df %>% filter(str_detect(Region,"Dist")) %>% mutate(Phenotype = fct_reorder(Phenotype, value) %>% fct_rev()) %>% pull(Phenotype) %>% levels()  

df %>% 
  mutate(Phenotype = fct_relevel(Phenotype, orderedLevels)) %>% 
  ggplot(aes(x= Region, y = value, alluvium = Phenotype))+
  geom_alluvium(aes(fill = Phenotype), alpha =0.75, decreasing = F, width = 0.02)+
  theme_void()+
  theme(axis.text.x = element_text(),
        plot.margin = margin(b = 10))

df <- combined %>% 
  select(PatientID, Region,contains("percent"),
         -matches("n_|Uni")) %>% 
  filter(!str_detect(Region,"total")) %>% 
  rename_all(~str_remove(.,"_percent")) %>% 
  pivot_longer(cols = where(is.numeric), names_to = c("Phenotype")) %>%
  group_by(Region, Phenotype) %>% 
  summarise(across(where(is.numeric),mean)) %>%
  filter(str_detect(Phenotype,"FAP|SMA|^FSP$")) %>% 
  ungroup() %>% 
  mutate(Region = str_to_title(Region) %>% str_replace(" ","\n"),
         Region = fct_relevel(Region,"Tumour", "Intratumoural\nStroma","Peritumoural\nStroma","Distal\nStroma"),
         Phenotype = fct_relevel(Phenotype),
         Phenotype = fct_reorder(Phenotype, value) %>% fct_rev())

orderedLevels <- df %>% filter(str_detect(Region,"Dist")) %>% mutate(Phenotype = fct_reorder(Phenotype, value) %>% fct_rev()) %>% pull(Phenotype) %>% levels()  

df %>% 
  mutate(Phenotype = fct_relevel(Phenotype, orderedLevels)) %>% 
  ggplot(aes(x= Region, y = value, alluvium = Phenotype))+
  geom_alluvium(aes(fill = Phenotype), alpha =0.75, decreasing = F, width = 0.02)+
  theme_void()+
  theme(axis.text.x = element_text(),
        plot.margin = margin(b = 10))

library(spatstat.geom)

getDistanceData <-function(data, Phenotype1,Phenotype2, k=1){
  idColString = ifelse(k==1, "which", paste0("which.",k))
  window <- owin(
    xrange = c(min(data$x), max(data$x)),
    yrange = c(min(data$y), max(data$y))
  )
  
  cellType1 <- data %>%
    filter(Phenotype == Phenotype1) %>%
    drop_na() 
  
  cellType2 <- data %>%
    filter(Phenotype == Phenotype2) %>%
    drop_na()
  
  distData <- nncross(ppp(cellType1$x,cellType1$y, window),
                      ppp(cellType2$x, cellType2$y, window),
                      k =k) %>% 
    rename_with(~str_remove(.,paste0("\\.",k))) %>% 
    mutate(target = Phenotype2)
  
  distData<-distData %>% 
    mutate(from.x = cellType1$x,
           from.y = cellType1$y,
           to.x = cellType2$x[which],
           to.y = cellType2$y[which]) %>% 
    cbind(cellType1 %>% select(Region,  PatientID, subtype, sampleType = `sample type`))
  
  
  distData
}
#allCells <- read_rds("Data/allCells.rds")
nestedCells <- allCells %>%
  filter(`sample type` == "resection") %>% 
  mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
  nest(.by = FOV)

# 
# distanceToCK<- allCells %>% 
#   mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
# getDistanceData("FAP","CK")
# 
# distanceToCD3<- allCells %>% 
#   mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
# getDistanceData("FAP","CD3")
# 
# distanceToCD14<- allCells %>% 
#   mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
# getDistanceData("FAP","CD14")
# 
# distanceToCD68<- allCells %>% 
#   mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
# getDistanceData("FAP","CD68")
# 
# distanceToCD68CD14<- allCells %>% 
#   mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
# getDistanceData("FAP","CD68+CD14")
# 
# distanceToCD68FSP<- allCells %>% 
#   mutate(Phenotype = ifelse(str_detect(Phenotype, "FAP"),"FAP", Phenotype)) %>% 
# getDistanceData("FAP","CD68+FSP")
# 
# distanceToCK %>%
#   ggplot(aes(x = dist)) +
#   geom_density(color = "red")+
#   geom_density(data = distanceToCD3, color = "blue")+
#   geom_density(data = distanceToCD14, color = "green")+
#   geom_density(data = distanceToCD68, color = "black")+
#   geom_density(data = distanceToCD68CD14, color = "magenta")+
#   geom_density(data = distanceToCD68FSP, color = "pink")+
#   scale_x_log10()+
#   theme_minimal()

# 
# distanceToCK<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CK") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCD3<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CD3") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCD14<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CD14") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCD68<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CD68") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCD68CD14<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CD68+CD14") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCD68FSP<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CD68+FSP") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCD68FSPCD14<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","CD68+FSP+CD14") %>% mutate(FOV = FOV)
#   })
# 
# distanceToFSPCD14<- nestedCells %>% 
#   pmap_df(function(FOV,data){
#     data %>% 
# getDistanceData("FAP","FSP+CD14") %>% mutate(FOV = FOV)
#   })
# 
# distanceToCK %>%
#   ggplot(aes(x = dist)) +
#   geom_density(color = "red")+
#   geom_density(data = distanceToCD3, color = "blue")+
#   geom_density(data = distanceToCD14, color = "green")+
#   geom_density(data = distanceToCD68, color = "black")+
#   geom_density(data = distanceToCD68CD14, color = "magenta")+
#   geom_density(data = distanceToCD68FSP, color = "pink")+
#   geom_density(data = distanceToCD68FSPCD14, color = "orange")+
#    geom_density(data = distanceToFSPCD14, color = "cyan")+
#   scale_x_log10()+
#   theme_minimal()

phenotypes = c("CK","CD68","CD14","CD3","CD68+CD14","CD68+FSP","CD68+FSP+CD14","FSP+CD14")

c <-distanceDataframe<- nestedCells %>% 
  pmap_df(function(FOV,data){
    
    map_df(phenotypes, function(phenotype)  {data %>% getDistanceData("FAP",phenotype) %>% mutate(FOV = FOV)})
  })

c %>%
  ggplot(aes(x = dist, color = target)) +
  geom_density()+
  scale_x_log10()+
  theme_minimal()


