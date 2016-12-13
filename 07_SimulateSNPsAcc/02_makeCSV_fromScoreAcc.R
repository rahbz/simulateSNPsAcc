args <- commandArgs(TRUE)
# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1 Path to HDF5 file of whole imputed dataset
pathHDF5file <- args[1]
#2 Input scoreacc file directory
#workingDir <- args[1]
#setwd(workingDir)
#working directory should be the same

#2 Output CSV
outFile <- args[2]
confScore = 0.95

library("rhdf5")
the1001accessions <- h5read(pathHDF5file, name = "accessions")

list <- list.files("./", pattern = "[.]ScoreAcc.txt$")
Names <- character()
ReqAcc <- character()
TopHitAcc <- character()
NextHitAcc <- character()
ScoresTop <- numeric()
FracEstimate <- numeric()
ScoreNOs <- numeric()
TopHitMatchedSNPs <- numeric()
NextHitMatchedSNPs <- numeric()
ChoiceAcc <- numeric()
pVal <- numeric()
SNPscalled <- numeric()
for (echscore in list){
  name <- sub("[.]ScoreAcc[.]txt", "", echscore)
  ScoreAcc <- read.table(echscore, header = F)
  givaccid <- as.numeric(unlist(strsplit(name, "_"))[1]) + 1
  snps <- as.numeric(unlist(strsplit(name, "_"))[2])
  givacc <- the1001accessions[givaccid]
  topone <- which(ScoreAcc$V4 == max(ScoreAcc$V4))
  nextacc <- ScoreAcc$V1[order(-ScoreAcc$V4)][2]
  nexthitmatchsnps <- ScoreAcc$V2[order(-ScoreAcc$V4)][2]
  #tune <- which(ScoreAcc$V4 > (max(ScoreAcc$V4)-2*sd(ScoreAcc$V4)))
  tune <- length(which(ScoreAcc$V5 > confScore))
  pVal <- c(pVal, ScoreAcc$V5[topone[1]])
  if(file.exists(paste(name, ".refScore.txt", sep = ""))){
    if(length(readLines(paste(name, ".refScore.txt", sep = ""))) > 0){
      refScore <- read.table(paste(name,".refScore.txt", sep = ""), header = F)
      tune <- length(which(refScore$V5 > confScore))
    }
  }
  
  choicenum <- which(ScoreAcc$V1[order(-ScoreAcc$V4)] == as.numeric(givacc))
  Names <- c(Names, name)
  ReqAcc <- c(ReqAcc, givacc)
#  print(topone)
  TopHitAcc <- c(TopHitAcc, ScoreAcc$V1[topone[1]])
  NextHitAcc <- c(NextHitAcc, nextacc)
  ScoresTop <- c(ScoresTop, max(ScoreAcc$V4))
  FracEstimate <- c(FracEstimate, ScoreAcc$V4[order(-ScoreAcc$V4)][2]/max(ScoreAcc$V4))
  ScoreNOs <- c(ScoreNOs, tune)
  TopHitMatchedSNPs <- c(TopHitMatchedSNPs, ScoreAcc$V2[topone[1]])
  NextHitMatchedSNPs <- c(NextHitMatchedSNPs, nexthitmatchsnps)
  ChoiceAcc <- c(ChoiceAcc, choicenum)
  SNPscalled <- c(SNPscalled, snps)
}

print(length(Names))
print(length(ReqAcc))
print(length(TopHitAcc))
print(length(NextHitAcc))
print(length(ScoresTop))
print(length(FracEstimate))
print(length(ScoreNOs))
print(length(pVal))
print(length(TopHitMatchedSNPs))
print(length(NextHitMatchedSNPs))
print(length(ChoiceAcc))
print(length(SNPscalled))

DF <- data.frame(FILENAME = Names, RequiredAccession = ReqAcc, TopHitAccession = TopHitAcc, NextHitAcc = NextHitAcc, Score = ScoresTop, FracEstimate = FracEstimate, MatchedSNPs = TopHitMatchedSNPs, NextHitMatchedSNPs = NextHitMatchedSNPs, SNPscalled = SNPscalled, pVal = pVal, TopHitsNumber = ScoreNOs, ChoiceNum = ChoiceAcc)

write.csv(DF, file = outFile)
