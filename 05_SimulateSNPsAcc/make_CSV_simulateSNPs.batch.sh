#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N makeCSV
#PBS -J 0-19
#PBS -V
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -o log.makeCSV
#PBS -e err.makeCSV


cd $PBS_O_WORKDIR

alldirs=(`ls -d */`)
hdf5=`ls wholeImputed_filtered_all_chroms.hdf5`

workDir="${alldirs[$PBS_ARRAY_INDEX]}"

module load R

Rscript ~/MyScripts/04_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R $PBS_O_WORKDIR/$workDir $PBS_O_WORKDIR/$workDir/intermediate_modified.csv $PBS_O_WORKDIR/$hdf5


#workDir=`pwd`
#
#for (( i=0; i < ${#alldirs[@]}; i++ ))
#do
#	Rscript ~/MyScripts/04_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R $workDir/${alldirs[$i]} $workDir/${alldirs[$i]}/intermediate_modified.csv $workDir/$hdf5 
##	if [ $(($i % 16)) = 0 ]
##	then
##		wait
##	fi
#
#done
#
#
