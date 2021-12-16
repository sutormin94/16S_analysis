#!bin/bash

PWD='/home/niagara/Storage/MetaRus/V_Mamontov/results_16S'
echo 'Start quality control!'

#Path to the file containing sequencing adapters sequences for trimmomatic uses. Typically in the Trimmomatic-0.36/adapters/XXX.fa
Adapters='/home/niagara/Progs/Trimmomatic-0.39/adapters/All_TruSeq.fa'
trimmomatic='/home/niagara/Progs/Trimmomatic-0.39/trimmomatic-0.39.jar'
samples='/home/niagara/Storage/MetaRus/V_Mamontov/16S_test_data/all_data'


#######
#Quality control and sequencing data preparation.
#######
echo '
#######################
Initial quality control is in progress...
#######################
'

#Initial quality control
mkdir $PWD/Fastqc_analysis/
mkdir $PWD/Fastqc_analysis/Initial
fastqc -t 40 -o $PWD/Fastqc_analysis/Initial $samples/*_R1.fastq.gz

#######
#Reads trimming
#######
echo '
#######################
Reads trimming...
#######################
'

mkdir $PWD/Trimmed
#for i in `ls -a $samples/ | grep 'fastq' | sed -r "s/(.+)_R[1,2]_001\.fastq\.gz/\1/g" | uniq | sort -d`; do
#for i in `ls -a $samples/ | grep '_F.fastq' | sed -r "s/(.+)_F\.fastq\.gz/\1/g" | uniq | sort -d`; do
for i in `ls -a $samples/ | grep '_R1.fastq' | sed -r "s/(.+)_R1\.fastq\.gz/\1/g" | uniq | sort -d`; do
echo $i
java -jar $trimmomatic SE -threads 40 -phred33 $samples/${i}_R1.fastq.gz $PWD/Trimmed/${i}_R1.fastq.gz ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 HEADCROP:17 MINLEN:150 ; done



#######
#Quality control after the trimming procedure
#######
echo '
#######################
Quality control after trimming...
#######################
'

mkdir $PWD/Fastqc_analysis/Trimmed/
fastqc -t 40 -o $PWD/Fastqc_analysis/Trimmed/ $PWD/Trimmed/*

cd $PWD/Trimmed

rm *unpaired_R2* | rm *R2* | rm *unpaired*

cd $PWD

