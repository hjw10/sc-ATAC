bash download.sh

prefetch -O ./1.data/ --option-file SRR_Acc_List.txt

fastq-dump --gzip --split-files ../1.data/SRR*/*

#改名字
#SRR18505563_S1_L001_R1_001.fastq.gz
#SRR18505563_S1_L001_R2_001.fastq.gz
#SRR18505563_S1_L001_R3_001.fastq.gz
#SRR18505564_S1_L001_R1_001.fastq.gz
#SRR18505564_S1_L001_R2_001.fastq.gz
#SRR18505564_S1_L001_R3_001.fastq.gz
#SRR18505565_S1_L001_R1_001.fastq.gz
#SRR18505565_S1_L001_R2_001.fastq.gz
#SRR18505565_S1_L001_R3_001.fastq.gz

#cellranger-atac
bin=~/tools/cellranger-atac-2.1.0/bin/cellranger-atac
ref=~/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0
for id in d30 d15 d8
do
time ${bin} count   --id ${id} \
                        --reference ${ref} \
                        --fastqs ./2.raw_fq/${id} \
                        --sample ${id} \
                        --localcores 8 \
                        --localmem 64
done

