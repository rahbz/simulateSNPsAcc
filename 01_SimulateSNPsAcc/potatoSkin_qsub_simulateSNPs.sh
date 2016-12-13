#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N simSNPs900
#PBS -V
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -o log

cd /lustre/scratch/users/rahul.pisupati/simulateSNPs/0.9kSNPs/
numSNPs=900
numAcc=2029
hdf5=/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5
hdf5_acc=/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5
totSNPs=/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt

module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

for (( i=0; i < $numAcc; i++ ))
do
	python ~/MyScripts/01_SimulateSNPsAcc/02_SimulateSNPsperAcc.py -d $hdf5 -e $hdf5_acc -t $totSNPs -a $i -n $numSNPs -o ${i}_$numSNPs.ScoreAcc.txt &
	if [ $(($i % 16)) = 0 ]
        then
                wait
        fi
done
wait
# Make a CSV file
module load R
Rscript ~/MyScripts/01_SimulateSNPsAcc/makeCSV_fromScoreAcc.R ./ intermediate_modified.csv

