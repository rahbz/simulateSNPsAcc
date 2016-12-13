#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N getPOSfile
#PBS -J 1-1135
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -o log.posVCF
#PBS -e err.posVCF

### The job array ID has to be 
### 2 -- (number_of_columns_in_CSV - 1)

cd $PBS_O_WORKDIR

CSVfile="/lustre/scratch/projects/the1001genomes/rahul/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.csv"
#CSVfile="/lustre/scratch/projects/the1001genomes/VCF_1163g_Jan20/1163g.179kB.prior15.gauss4.ts99.5.BIALLELIC.BINARY.csv"
#CSVfile="/lustre/scratch/projects/the1001genomes/VCF_1163g/Swedes243.177k.prior15.gauss6.ts99.5.BIALLELIC.hetfiltered.csv"
#CSVfile="/lustre/scratch/projects/field_experiments/007.pilot.sequencing/018.genotyping.by.plate/997.Swedes.220.10May2016/02_2.3M_Swedes220.175k.prior15.gauss4.ts99.5.BIALLELIC.csv"
#CSVfile="/lustre/scratch/projects/aquilegia/998.rahul_snpmatch/all.49inds.plus.semiaq.filter.multi.indel.repeat.BIALLELIC.csv"

#outFol="posFiles"
#mkdir $outFol

#list=(`cut -d "," -f 1 /lustre/scratch/projects/the1001genomes/rahul/check_Tnz-0/A_thaliana_sequencing_whishlist-filtered.csv |tail -n +2`)
#list=(`cat check.txt`)
#id=`head -n 1 $CSVfile | sed 's/,/\n/g' | grep -nw "${list[$PBS_ARRAY_INDEX]}" | cut -f1 -d":"`
#awk -F"," -v nu=$id '$nu == 0 {print $1 "\t" $2 "\t0/0"} $nu == 1 {print $1 "\t" $2 "\t1/1"}'  $CSVfile > ${list[$PBS_ARRAY_INDEX]}.pos.txt 

## get POS files for all the accessions
list=(`head -n 1 $CSVfile | sed 's/,/\n/g'`)
#awk -F"," -v nu=$((PBS_ARRAY_INDEX + 2)) '$nu == 0 {print $1 "\t" $2 "\t0/0"} $nu == 1 {print $1 "\t" $2 "\t1/1"} $nu == 2 {print $1 "\t" $2 "\t0/1"}' $CSVfile > $outFol/${list[$((PBS_ARRAY_INDEX + 1))]}.pos.txt 
awk -F"," -v nu=$((PBS_ARRAY_INDEX + 2)) '$nu == 0 {print $1 "\t" $2 "\t0/0"} $nu == 1 {print $1 "\t" $2 "\t1/1"} $nu == 2 {print $1 "\t" $2 "\t0/1"}' $CSVfile > ${list[$((PBS_ARRAY_INDEX + 1))]}.pos.txt 

#for i in `cut -d "," -f 1 /lustre/scratch/projects/the1001genomes/rahul/check_Tnz-0/A_thaliana_sequencing_whishlist-filtered.csv |tail -n +2`
#do
#	id=`head -n 1 $CSVfile | sed 's/,/\n/g' | grep -nw "$i" | cut -f1 -d":"`
#	awk -F"," -v nu=$id '$nu == 0 {print $1 "\t" $2 "\t0/0"} $nu == 1 {print $1 "\t" $2 "\t1/1"}'  $CSVfile > $i.pos.txt
#done







