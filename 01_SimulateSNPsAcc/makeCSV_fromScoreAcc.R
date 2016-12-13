args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1 Path to HDF5 file of whole imputed dataset
pathHDF5file <- "/lustre/scratch/users/rahul.pisupati/wholeImputed_all_chromosomes_binary_acc.hdf5"
#2 Input scoreacc file directory
workingDir <- args[1]
setwd(workingDir)

#2 Output CSV
outFile <- args[2]


library("rhdf5")
the1001accessions <- h5read(pathHDF5file, name = "accessions")

list <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
ReqAcc <- character()
TopHitAcc <- character()
ScoresTop <- character()
ScoreNOs <- character()
ScoreSDs <- character()
TopHitMatchedSNPs <- character()
ChoiceAcc <- character()
for (echscore in list){
  name <- sub("[.]ScoreAcc[.]txt", "", echscore)
  ScoreAcc <- read.table(echscore, header = F)
  givaccid <- as.numeric(unlist(strsplit(name, "_"))[1]) + 1
  givacc <- the1001accessions[givaccid]
  topone <- which(ScoreAcc$V4 == max(ScoreAcc$V4))
  tune <- which(ScoreAcc$V4 > (max(ScoreAcc$V4)-2*sd(ScoreAcc$V4)))
  tstat <- (max(ScoreAcc$V4) - mean(ScoreAcc$V4))/sd(ScoreAcc$V4)
  choicenum <- which(ScoreAcc$V1[order(-ScoreAcc$V4)] == as.numeric(givacc))
  Names <- c(Names, name)
  ReqAcc <- c(ReqAcc, givacc)
  TopHitAcc <- c(TopHitAcc, ScoreAcc$V1[topone])
  ScoresTop <- c(ScoresTop, max(ScoreAcc$V4))
  ScoreNOs <- c(ScoreNOs, length(tune))
  ScoreSDs <- c(ScoreSDs, tstat)
  TopHitMatchedSNPs <- c(TopHitMatchedSNPs, ScoreAcc$V2[topone])
  ChoiceAcc <- c(ChoiceAcc, choicenum)
}
DF <- data.frame(FILENAME = Names, RequiredAccession = ReqAcc, TopHitAccession = TopHitAcc, Score = as.numeric(ScoresTop), MatchedSNPs = as.numeric(TopHitMatchedSNPs), tStat  = as.numeric(ScoreSDs), TopHitsNumber = ScoreNOs, ChoiceNum = ChoiceAcc)

write.csv(DF, file = outFile)

