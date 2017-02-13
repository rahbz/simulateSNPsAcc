#!/bin/bash
#PBS -S /bin/bash
#PBS -P the1001genomes
#PBS -N sim_SNPmatch
#PBS -q new_nodes
#PBS -J 1-985
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=8:mem=80gb
#PBS -o log.simSNPs
#PBS -e err.simSNPs

cd $PBS_O_WORKDIR
nct=8

#hdf5="/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.hdf5"
#hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/1135g_SNP_BIALLELIC.hetfiltered.snpmat.6oct2015.acc.hdf5"
#accList="/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/the1001accessions.list.txt"

hdf5="/lustre/scratch/projects/the1001genomes/rahul/108.simulated.1001g.SNPmatrix/002.num_5k/the1001genomes_filtered_all_chroms.hdf5"
hdf5_acc="/lustre/scratch/projects/the1001genomes/rahul/108.simulated.1001g.SNPmatrix/002.num_5k/the1001genomes_filtered_all_chroms.acc.hdf5"
accList="/lustre/scratch/projects/the1001genomes/rahul/108.simulated.1001g.SNPmatrix/002.num_5k/the1001genomes_filtered_all_chroms.accList.txt"

# wholeImputed SNP matrix
#hdf5=`readlink -f wholeImputed_filtered_all_chroms.hdf5`
#hdf5_acc=`readlink -f wholeImputed_filtered_all_chroms.acc.hdf5`

#  SWEDES SNP matrix
#hdf5=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.hdf5
#hdf5_acc=/lustre/scratch/projects/the1001genomes/rahul/SWEDEs/Swede.198.hetfiltered.bialleleic.22Jun15.binary.hdf5

module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
module load snpmatch


accID=`cat $accList | head -n $PBS_ARRAY_INDEX | tail -n 1`

nSNPs=(`echo "100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000,6000000,7000000,8000000,9000000" | sed 's/,/\n/g'`)
numSims=${#nSNPs[@]}


for (( i=0; i < $numSims; i++ ))
do
	mkdir -p SNPs_${nSNPs[$i]}
	mkdir -p logs/logs_SNPs_${nSNPs[$i]}
	python ~/MyScripts/simulateSNPsAcc/08_SimulateSNPsAcc/01_SimulateSNPsperAcc.py  -d $hdf5 -e $hdf5_acc -a $accID -n ${nSNPs[$i]} -o SNPs_${nSNPs[$i]}/${accID}_${nSNPs[$i]}.ScoreAcc.txt > logs/logs_SNPs_${nSNPs[$i]}/logs.${accID} 2>&1 &
	if [ $(($i % 8)) == 0 ]
	then
		wait
	fi
done
wait
