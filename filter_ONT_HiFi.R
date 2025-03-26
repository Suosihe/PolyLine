library(tidyverse)

setwd("D:/R")

ont <- read_tsv("all_ONT.tsv", col_names = F)

ont_select <- ont %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X4)), x2 = min_rank(desc(X3))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 %in% c("HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", 
                   "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6", 
                   "HiC_scaffold_7", "HiC_scaffold_8", "HiC_scaffold_9", 
                   "HiC_scaffold_10", "HiC_scaffold_12", "HiC_scaffold_14", 
                   "HiC_scaffold_15", "HiC_scaffold_16", "HiC_scaffold_18", 
                   "HiC_scaffold_19", "HiC_scaffold_20")) %>% 
  filter(X4 != 0)
  

write_tsv(ont_select, "select.ONT.tsv", col_names = T)

hifi <- read_tsv("all_HiFi.tsv", col_names = F)

hifi_select <- hifi %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X4)), x2 = min_rank(desc(X3))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 %in% c("HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", 
                   "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6", 
                   "HiC_scaffold_7", "HiC_scaffold_8", "HiC_scaffold_9", 
                   "HiC_scaffold_10", "HiC_scaffold_12", "HiC_scaffold_14", 
                   "HiC_scaffold_15", "HiC_scaffold_16", "HiC_scaffold_18", 
                   "HiC_scaffold_19", "HiC_scaffold_20")) %>% 
  filter(X4 != 0)

write_tsv(hifi_select, "select.HiFi.tsv", col_names = T)
