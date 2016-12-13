#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N simSNPs
#PBS -J 0-21
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=42gb
#PBS -o log.simSNPs
#PBS -e err.simSNPs

cd $PBS_O_WORKDIR

#hdf5="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
#kinFile="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.kinship.ibs.hdf5"
#totSNPs="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.numInformative.txt"
#totSNPs="/lustre/scratch/projects/the1001genomes/VCF_1135g/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.totSNP.count.txt"

hdf5=`ls the1001genomes_filtered_all_chroms.hdf5`
hdf5_acc=`ls the1001genomes_filtered_all_chroms.acc.hdf5`
kinFile=`ls the1001genomes_filtered_all_chroms.kinship.ibs.hdf5`
totSNPs=`ls the1001genomes_filtered_all_chroms.numInformative.txt`
#totSNPs=`ls the1001genomes_filtered_all_chroms.totalSNPcount.txt`

# wholeImputed SNP matrix
#hdf5=`ls wholeImputed_filtered_all_chroms.hdf5`
#hdf5_acc=`ls wholeImputed_filtered_all_chroms.acc.hdf5`
#kinFile=`ls wholeImputed_filtered_all_chroms.kinship.ibs.hdf5`
#totSNPs=`ls wholeImputed_filtered_all_chroms.totalSNPcount.txt`

#  SWEDES SNP matrix
#hdf5=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.hdf5
#hdf5_acc=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.hdf5
#totSNPs=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.totalAccSNPcount.txt
#kinFile=/lustre/scratch/users/rahul.pisupati/wholeImputed_kinship_ibs_binary_mac5.h5py

module load h5py
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

numAcc=`wc -l $totSNPs | cut -f1 -d ' '`
nSNPs=(`echo "100,200,300,400,500,800,1000,2000,3000,4000,5000,8000,10000,20000,30000,50000,80000,100000,200000,300000,400000,500000" | sed 's/,/\n/g'`)

nSNP=${nSNPs[$PBS_ARRAY_INDEX]}
mkdir SNPs_$nSNP
cd SNPs_$nSNP

for (( i=0; i < $numAcc; i++ ))
do
	python ~/MyScripts/simulateSNPsAcc/05_SimulateSNPsAcc/01_SimulateSNPsperAcc.py -d $PBS_O_WORKDIR/$hdf5 -e $PBS_O_WORKDIR/$hdf5_acc -t $PBS_O_WORKDIR/$totSNPs -a $i -n $nSNP -o ${i}_$nSNP.ScoreAcc.txt -r ${i}_$nSNP.refScore.txt &
#	python ~/MyScripts/simulateSNPsAcc/05_SimulateSNPsAcc/01_SimulateSNPsperAcc.py -d $hdf5 -e $hdf5_acc -t $totSNPs -a $i -n $nSNP -o ${i}_$nSNP.ScoreAcc.txt -r ${i}_$nSNP.refScore.txt &
	if [ $(($i % 16)) = 0 ]
	then
		wait
	fi
done
wait

# Make a CSV file
module load R
#Rscript ~/MyScripts/simulateSNPsAcc/05_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R $PBS_O_WORKDIR/SNPs_$nSNP $PBS_O_WORKDIR/SNPs_$nSNP/intermediate_modified.csv $hdf5
Rscript ~/MyScripts/simulateSNPsAcc/05_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R $PBS_O_WORKDIR/SNPs_$nSNP $PBS_O_WORKDIR/SNPs_$nSNP/intermediate_modified.csv $PBS_O_WORKDIR/$hdf5
