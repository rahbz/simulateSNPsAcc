#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N simSNPs500k
#PBS -V
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -o log


cd $PBS_O_WORKDIR
numSNPs=500000

hdf5=`ls ../the1001genomes_filtered_all_chroms.hdf5`
hdf5_acc=`ls ../the1001genomes_filtered_all_chroms.acc.hdf5`
kinFile=`ls ../the1001genomes_filtered_all_chroms.kinship.ibs.hdf5`
totSNPs=`ls ../the1001genomes_filtered_all_chroms.totalSNPcount.txt`

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.kinship.ibs.hdf5"
#totSNPs="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.totSNP.count.txt"

#hdf5=/lustre/scratch/users/rahul.pisupati/the1001genomes_all_chromosomes_binary.hdf5
#hdf5_acc=/lustre/scratch/users/rahul.pisupati/the1001genomes_all_chromosomes_binary_acc.hdf5
#totSNPs=/lustre/scratch/users/rahul.pisupati/the1001genomes_totalSNPs_peracc.txt

numAcc=`wc -l $totSNPs | cut -f1 -d ' '`
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

for (( i=0; i < $numAcc; i++ ))
do
	python ~/MyScripts/simulateSNPsAcc/04_SimulateSNPsAcc/01_SimulateSNPsperAcc.py -d $hdf5 -e $hdf5_acc -t $totSNPs -a $i -n $numSNPs -o ${i}_$numSNPs.ScoreAcc.txt -r ${i}_$numSNPs.refScore.txt &
	if [ $(($i % 16)) = 0 ]
        then
                wait
        fi
done
wait
# Make a CSV file
module load R
#Rscript ~/MyScripts/02_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R ./ intermediate_modified.csv $numSNPs
Rscript ~/MyScripts/simulateSNPsAcc/04_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R ./ intermediate_modified.csv $hdf5
