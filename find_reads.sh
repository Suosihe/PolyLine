#!/bin/bash

## `euulhl0.hic.p_ctg.FINAL.fasta` is the preliminary result obtained after hifiasm and juicer pipeline

## for ONT and HiFi reads

minimap2 -x map-ont euulhl0.hic.p_ctg.FINAL.fasta SK-1-ONT/pass.all.fq.gz -t 40 > all_ONT.paf

minimap2 -x map-hifi euulhl0.hic.p_ctg.FINAL.fasta SK-1-HiFi/hifi_reads/sk_hifi.fasta -t 40 > all_HiFi.paf

cut -f 1,6,10,12 all_ONT.paf>all_ONT.tsv
cut -f 1,6,10,12 all_HiFi.paf>all_HiFi.tsv

# cut -f 1,6,8,9,10,12 all_ONT.paf > all_ONT2.tsv
# cut -f 1,6,8,9,10,12 all_HiFi.paf > all_HiFi2.tsv

Rscript filter_ONT_HiFi.R

for i in {1..20};do echo "awk '{if($2==\"HiC_scaffold_$i\"){print}}' select.ONT.tsv |cut -f 1 > ONT_$i.name";done >> getfastq.sh
for i in {1..20};do echo "awk '{if($2==\"HiC_scaffold_$i\"){print}}' select.HiFi.tsv |cut -f 1 > HiFi_$i.name";done >> getfastq.sh


## for Hi-C reads

cut -d " " -f 2,15 merged_nodups.txt |sed 's/ /\t/g' > ctg-hic.tsv

### 該腳本需要修訂：1-2的內容 如何設置參數？
Rscript 3-17names.R

for i in *hic.name; do gawk '{if ($2 != "X2") {print $2 "/1"}}' $i > $i.R1; done
for i in *hic.name; do gawk '{if ($2 != "X2") {print $2 "/2"}}' $i > $i.R2; done

for i in *hic.name.R1; do echo "seqtk subseq SK-hic_R1.fastq $i > $i.fastq"; done >> getfastq2.sh
for i in *hic.name.R2; do echo "seqtk subseq SK-hic_R2.fastq $i > $i.fastq"; done >> getfastq2.sh

## select the Scaffold specific reads

sh getfastq.sh
sh getfastq2.sh
