#!/usr/bin/python
import math
import random
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs to be taken for analysis", type="int")
inOptions.add_option("-a", "--acc_to_check", dest="acc", help="Index of the accession needed to be checked", type="int")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")

(options, args) = inOptions.parse_args()

GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

AccToCheck = options.acc
numSNPs = options.numsnps

sampleSNPs = numpy.sort(random.sample(numpy.where(GenotypeData_acc.snps[:,AccToCheck] == 1)[0], numSNPs))
ScoreList = numpy.zeros(len(GenotypeData.accessions))
for i in range(0, len(sampleSNPs), 1000):
  subsampleSNPs = tuple(sampleSNPs[i:i+1000])
  ScoreList = ScoreList + numpy.sum(GenotypeData.snps[subsampleSNPs,], axis = 0)

# Calculate the number of SNPs in all the accessions
filetotal = open(options.file_num_snps, 'r')
TotnumberSNPs = filetotal.read().split("\n")
filetotal.close()

#FinalScore = numpy.zeros(len(GenotypeData.accessions))
#FinalScore = [ScoreList[i]/float(TotnumberSNPs[i].split("\t")[1]) for i in range(0, len(GenotypeData.accessions))]

outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = TotnumberSNPs[i].split("\t")[1]
  outScore = (float(ScoreList[i])*float(ScoreList[i]))/(int(numsnp)*numSNPs)
  outfile.write(GenotypeData.accessions[i])
  outfile.write("\t")
  outfile.write("%s" % int(ScoreList[i]))
  outfile.write("\t")
  outfile.write("%s" % int(numsnp))
  outfile.write("\t")
  outfile.write("%s\n" % outScore)
outfile.close()



#TopHit = numpy.where(FinalScore == max(FinalScore))[0][0]
#TopHitsNum = len(numpy.where(FinalScore > (max(FinalScore) - (2*numpy.std(FinalScore))))[0])
#tStat = (max(FinalScore) - numpy.mean(FinalScore))/numpy.std(FinalScore)
#if TopHitsNum == 2:
#  print TopHit, "\t", TopHitsNum, "\t", max(FinalScore), "\t",  tStat
