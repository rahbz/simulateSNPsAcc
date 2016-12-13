# Generate a CSV file based on the ScoreAcc files generated before
# Input: 
# Various variables:
#1) Path to HDF5 file of whole imputed dataset
#pathHDF5file <- args[1]
#2) Run the script from the same working directory as the ScoreAcc files
#workingDir <- "./"
#setwd(workingDir)
#4) Path to the accession file, where the accessions supposed to be in a tab-delimited file
#5) Path to merged ecotype IDs
#6) Output CSV file

library("optparse")
option_list = list(
  make_option(c("-o", "--outFile"), type="character", default="intermediate_modified.csv", help="Output file name [default= %default]", metavar="character"),
  make_option(c("-f", "--folID"), type="character", default=NULL, help="Plate ID to be added in the CSV file", metavar="character")
);
opt = parse_args(OptionParser(option_list=option_list));


outFile = opt$outFile
folid = opt$folID

#library("rhdf5")
#the1001positions <- h5read(pathHDF5file, name="positions", read.attributes=T);
#the1001chromosomes <- attributes(the1001positions)$chr_regions;
#the1001chromosomes[1,] <- the1001chromosomes[1,] + 1;
#the1001chr.ids <- attributes(the1001positions)$chrs;
#the1001accessions <- h5read(pathHDF5file, name="accessions")

if(file.exists(outFile)){
  file.remove(outFile)
}

#-------
allScoreFiles <- list.files("./", pattern = "[.]ScoreAcc.txt$")
AssignedAcc <- character()
TopHitAcc <- character()
NextHitAcc <- character()
ThirdHit <- character()
TopHitScore <- numeric()
FracScore <- numeric()
LikeLihoodTopHit <- numeric()
LLRNextHit <- numeric()
TopHitsNumber <- numeric()
TopHits <- character()
ChoiceAcc <- numeric()
for (file in allScoreFiles){
  ScoreAcc <- try(read.table(file, header = F))
  if(inherits(ScoreAcc, 'try-error')){
    next
  }

  name <- sub(".ScoreAcc.txt","",file)
  folid <- unlist(strsplit(name, "_"))[2]
  assacc <- unlist(strsplit(name, "_"))[1]
  AssignedAcc <- c(AssignedAcc, assacc)
## Changed the sorting order based on the score to the likelihood ratio
#  ranks <- order(-ScoreAcc$V4, na.last = TRUE, partial = ScoreAcc$V6)
  ranks <- order(ScoreAcc$V6, na.last = TRUE, partial = ScoreAcc$V4)

  topscore <- ScoreAcc$V4[ranks[1]]
  topacc <- as.character(ScoreAcc$V1[ranks[1]])
  maxlike <- ScoreAcc$V5[ranks[1]]

  newLike <- ScoreAcc$V5[ranks]/ScoreAcc$V5[ranks[1]]
  nextlike <- newLike[2]
  
#  snps <- as.numeric(numSNPs$V1[which(numSNPs$V2 == name)])
  nextacc <- as.character(ScoreAcc$V1[ranks[2]])
  nextscore <- ScoreAcc$V4[ranks[2]]
  frac <- nextscore/topscore

  thirdHit <- as.character(ScoreAcc$V1[ranks[3]])
  topnum <- length(which(ScoreAcc$V6 < 3.841))
  TopHitsNumber <- c(TopHitsNumber, topnum)
  if(topnum > 20){
    TopHits <- c(TopHits, "TooManytoPrint")
  } else if(topnum > 3) {
    alltops <- paste(ScoreAcc$V1[ranks[4:topnum]], collapse = ":")
    TopHits <- c(TopHits, alltops)
  } else {
    TopHits <- c(TopHits, "NA")
  }

  choiceacc <- which(ScoreAcc$V1[ranks] == assacc)
  ChoiceAcc <- c(ChoiceAcc, choiceacc)

  FracScore <- c(FracScore, frac)
  TopHitAcc <- c(TopHitAcc, topacc)
  NextHitAcc <- c(NextHitAcc, nextacc)
  ThirdHit <- c(ThirdHit, thirdHit)
  TopHitScore <- c(TopHitScore, topscore)
  LikeLihoodTopHit <- c(LikeLihoodTopHit, maxlike)
  LLRNextHit <- c(LLRNextHit, nextlike)
}

fol = rep(folid, length(AssignedAcc))

DF <- cbind(FOL = fol, AssignedAcc = AssignedAcc, TopHitAccession = TopHitAcc, NextHit = NextHitAcc, ThirdHit = ThirdHit, Score = as.numeric(TopHitScore), FracScore = as.numeric(FracScore), LikelihoodRatio = LikeLihoodTopHit, NextHitLLR = LLRNextHit, TopHitsNumber = TopHitsNumber, ChoiceAcc = ChoiceAcc, TopHits = TopHits)

#----

write.csv(DF, file = outFile, quote = F)



