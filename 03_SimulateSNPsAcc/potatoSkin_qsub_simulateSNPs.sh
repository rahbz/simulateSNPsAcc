#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N simSNPs
#PBS -J 0-1661
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=2:mem=5gb
#PBS -o log.simSNPs
#PBS -e err.simSNPs

cd $PBS_O_WORKDIR
numSNPs=500000
#numAcc=1662
#numAcc=198
#numAcc=2029

hdf5="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_all_chromosomes_filtered.hdf5"
hdf5_acc="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_all_chromosomes_filtered_acc.hdf5"
totSNPs="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_totalSNPs_peracc.txt"
kinFile="/lustre/scratch/users/rahul.pisupati/simulateSNPs/wholeImputed_kinship_filtered.hdf5"

#hdf5=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.hdf5
#hdf5_acc=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.hdf5
#totSNPs=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.totalAccSNPcount.txt
#kinFile=/lustre/scratch/users/rahul.pisupati/wholeImputed_kinship_ibs_binary_mac5.h5py

module load h5py
module load numpy
module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load pygwas
#rm avgKinScore_$numSNPs.txt

#python ~/MyScripts/03_SimulateSNPsAcc/01_SimulateSNPsperAcc.py -d $hdf5 -e $hdf5_acc -t $totSNPs -a $i -n $numSNPs -o ${i}_$numSNPs.ScoreAcc.txt -k $kinFile -s avgKinScore_$numSNPs.txt 
python ~/MyScripts/03_SimulateSNPsAcc/01_SimulateSNPsperAcc.py -d $hdf5 -e $hdf5_acc -t $totSNPs -a $PBS_ARRAY_INDEX -n $numSNPs -o ${PBS_ARRAY_INDEX}_$numSNPs.ScoreAcc.txt -r ${PBS_ARRAY_INDEX}_$numSNPs.refScoreAcc.txt 

# Make a CSV file
#module load R
#Rscript ~/MyScripts/03_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R ./ intermediate_modified.csv $hdf5
