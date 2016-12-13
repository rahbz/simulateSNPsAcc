#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N simSNPs500
#PBS -V
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=32gb
#PBS -o log

cd /lustre/scratch/users/rahul.pisupati/simulateSNPs/0.5kSNPs/
numSNPs=500
numAcc=2029


hdf5=/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary.hdf5
hdf5_acc=/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5
totSNPs=/lustre/scratch/users/rahul.pisupati/wholeImputed_totalSNPs_peracc.txt
kinFile=/lustre/scratch/users/rahul.pisupati/wholeImputed_kinship_ibs_binary_mac5.h5py

module load h5py
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas

for (( i=0; i<$numAcc; i++ ))
do
	python ~/MyScripts/02_SimulateSNPsAcc/03_GetKinScoreAmbiguousAcc.py -n $numSNPs -a $i -k $kinFile -s ${i}_$numSNPs.ScoreAcc.txt -o avgKinScore_$numSNPs.txt &
	if [ $(($i % 16)) = 0 ]
	then
		wait
	fi
done
wait



