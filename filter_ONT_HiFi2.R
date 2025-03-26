library(tidyverse)

setwd("D:/R")

rm(list=ls())

ont <- read_tsv("all_ONT.tsv", col_names = F)

ont2 <- read_tsv("all_ONT2.tsv", col_names = F)

ont_sel2 <- ont %>% filter(X4 >= 48) %>% 
  filter(X2 == "HiC_scaffold_1" | X2 == "HiC_scaffold_2") %>% 
  select(X1)
  
ont_sel3 <- ont %>% filter(X4 >= 48) %>% 
  filter(X2 == "HiC_scaffold_15") %>% 
  select(X1)

write_tsv(ont_sel2, "ONT_1-2.name3", col_names = F)

write_tsv(ont_sel3, "ONT_15.name3", col_names = F)

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
  
ont_select2 <- ont2 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 %in% c("HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", 
                   "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6", 
                   "HiC_scaffold_7", "HiC_scaffold_8", "HiC_scaffold_9", 
                   "HiC_scaffold_10", "HiC_scaffold_12", "HiC_scaffold_14", 
                   "HiC_scaffold_15", "HiC_scaffold_16", "HiC_scaffold_18", 
                   "HiC_scaffold_19", "HiC_scaffold_20")) %>% 
  filter(X6 != 0)

ont_count <- ont_select %>% 
  group_by(X2) %>% 
  summarise(n = sum(X3))

write_tsv(ont_count, "ONT_count.tsv", col_names = T)

write_tsv(ont_select2, "select2.ONT.tsv", col_names = T)

write_tsv(ont_select, "select.ONT.tsv", col_names = T)

hifi <- read_tsv("all_HiFi.tsv", col_names = F)

hifi2 <- read_tsv("all_HiFi2.tsv", col_names = F)

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

hifi_select2 <- hifi2 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 %in% c("HiC_scaffold_1", "HiC_scaffold_2", "HiC_scaffold_3", 
                   "HiC_scaffold_4", "HiC_scaffold_5", "HiC_scaffold_6", 
                   "HiC_scaffold_7", "HiC_scaffold_8", "HiC_scaffold_9", 
                   "HiC_scaffold_10", "HiC_scaffold_12", "HiC_scaffold_14", 
                   "HiC_scaffold_15", "HiC_scaffold_16", "HiC_scaffold_18", 
                   "HiC_scaffold_19", "HiC_scaffold_20")) %>% 
  filter(X6 != 0)

write_tsv(hifi_select2, "select2.HiFi.tsv", col_names = T)

hifi_sel2 <- hifi %>% filter(X4 >= 48) %>% 
  filter(X2 == "HiC_scaffold_1" | X2 == "HiC_scaffold_2") %>% 
  select(X1)

hifi_sel3 <- hifi %>% filter(X4 >= 48) %>% 
  filter(X2 == "HiC_scaffold_15") %>% 
  select(X1)

hifi_sel18 <- hifi2 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 == "HiC_scaffold_18") %>% 
  filter(X6 != 0) 

write_tsv(hifi_sel2, "HiFi_1-2.name3", col_names = F)

write_tsv(hifi_sel3, "HiFi_15.name3", col_names = F)

hifi_count <- hifi_select %>% 
  group_by(X2) %>% 
  summarise(n = sum(X3))

write_tsv(hifi_count, "HiFi_count.tsv", col_names = T)

# minimap2 -x map-hifi L0.review.genome.FINAL.fasta ce_hifi.fq -t 20 -o ce_hifi.paf
# cut -f 1,6,8,9,10,12 ce_hifi.paf> ce_HiFi.tsv

ce_hifi <- read_tsv("ce_HiFi.tsv", col_names = F)

ce_hifi_select2 <- ce_hifi %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X4)), x2 = min_rank(desc(X3))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 %in% c("HiC_scaffold_2")) %>% 
  filter(X6 != 0)

write_tsv(ce_hifi_select2, "ce_hifi2.name", col_names = F)

ce_ont <- read_tsv("ce_ont.tsv", col_names = F)

ce_ont_select2 <- ce_ont %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X4)), x2 = min_rank(desc(X3))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 %in% c("HiC_scaffold_2")) %>% 
  filter(X6 != 0)

write_tsv(ce_ont_select2, "ce_ont.name", col_names = F)

hifi3 <- read_tsv("all_HiFi3.tsv", col_names = F)
ont3 <- read_tsv("all_ONT3.tsv", col_names = F)

xq_hifi <- hifi3 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 == "Hs_18_1") %>% 
  filter(X3 < 4600000 & X4 < 4600000) %>% 
  filter(X6 > 10)

xq_hifi2 <- hifi3 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 == "Hs_18_2") %>% 
  filter(X3 < 6810000 & X4 < 6810000) %>% 
  filter(X6 > 10)

xq_ont <- ont3 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 == "Hs_18_1") %>% 
  filter(X3 < 4600000 & X4 < 4600000) %>% 
  filter(X6 > 10)

xq_ont2 <- ont3 %>% group_by(X1) %>% 
  mutate(x1 = min_rank(desc(X6)), x2 = min_rank(desc(X5))) %>% 
  filter(x1 == 1 & x2 == 1) %>% 
  filter(X2 == "Hs_18_2") %>% 
  filter(X3 < 6810000 & X4 < 6810000) %>% 
  filter(X6 > 10)

write_tsv(xq_hifi, "xq_hifi1.tsv", col_names = T)
write_tsv(xq_hifi2, "xq_hifi2.tsv", col_names = T)
write_tsv(xq_ont, "xq_ont1.tsv", col_names = T)
write_tsv(xq_ont2, "xq_ont2.tsv", col_names = T)
