
## Load modules R
#module load R

ecos="/lustre/scratch/projects/the1001genomes/rahul/101.VCF_1001G_1135/ecotypesids_merged_added.csv"

rootFol=`pwd`
folders=(`ls -d $rootFol/SNPs_* `)

script="~/MyScripts/simulateSNPsAcc/08_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R"

length=${#folders[@]}
for (( i=0; i<$length ;i=i+1 ));do
	cd ${folders[$i]}/
	echo ${folders[$i]}
	Rscript ~/MyScripts/simulateSNPsAcc/08_SimulateSNPsAcc/02_makeCSV_fromScoreAcc.R 
done

cd $rootFol

#(head -n 1 ${IntFiles[0]} && tail -n +2 ${IntFiles[@]} ) | grep -v "^=" | sed '/^\s*$/d' 

