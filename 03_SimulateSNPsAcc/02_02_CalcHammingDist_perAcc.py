#!/usr/bin/python
import math
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
from pygwas.core import genotype


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-r", "--req_acc_list", dest="req_acc_list", help="Input only the required accessions to the script", type="string")
inOptions.add_option("-p", "--vcf_pos", dest="posFile", help="Input Position File trimmed from the VCF based on quality", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-d", "--acc_hdf5_file", dest="hdf5File", help="HDF5 file chunked along the columns", type="string")

(options, args) = inOptions.parse_args()


GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

# Reading the required accession file into a tuple

reqacc = open(options.req_acc_list, 'r')
ReqAccs= tuple(reqacc.read().rstrip().split(","))
AccList = tuple(numpy.where(numpy.in1d(GenotypeData.accessions, ReqAccs))[0])

# Find all the positions where the AccList differs

AccSNPs = GenotypeData.snps[:,AccList]
max = len(AccList)
sumNumpy = numpy.sum(AccSNPs,axis = 1)
RefPosInd = numpy.where((sumNumpy != 0) & (sumNumpy != max))[0]

totalAccSNPs = numpy.sum(AccSNPs[RefPosInd,:], axis = 0)

# Now we have the SNPs where Accessions differ
# Take pos.txt file and check how many of them are matching with the SNPs 
# Take into accout the chromosome regions.
#
posfile = open(options.posFile, 'r')
TargetSNPs = tuple(posfile.read().rstrip().split("\n"))
ScoreList = numpy.zeros(len(AccList))
for i in range(0, len(GenotypeData.chr_regions)):
  start = GenotypeData.chr_regions[i][0]
  end = GenotypeData.chr_regions[i][1]
  chrNo = i + 1
  perChrPos = GenotypeData.positions[start:end]
  perChrTargetSNPs = [float(j.split("\t")[1]) for j in TargetSNPs if int(j.split("\t")[0]) == chrNo]
  matchedTargetSNPind = numpy.where(numpy.in1d(perChrPos, perChrTargetSNPs)) + start
# RefPosition Indices
  perChrRefPosInd = RefPosInd[numpy.where((RefPosInd > start) & (RefPosInd < end))]
# Take all the matched indices
  matRefPosInd = perChrRefPosInd[numpy.where(numpy.in1d(perChrRefPosInd, matchedTargetSNPind))]
# Now add the scores only for these matched Ref Pos Ind
  ScoreList = ScoreList + numpy.sum(AccSNPs[matRefPosInd,:], axis = 0)
  
# Number of SNPs that actually differ in the accessions

outfile = open(options.outFile, 'w')
for i in range(0, len(AccList)):
  numSNPAcc = totalAccSNPs[i]
  outScore = float(ScoreList[i])/numSNPAcc
  acc = GenotypeData.accessions[AccList[i]]
  outfile.write("%s" % int(acc))
  outfile.write("\t")
  outfile.write("%s" % int(ScoreList[i]))
  outfile.write("\t")
  outfile.write("%s" % int(numSNPAcc))
  outfile.write("\t")
  outfile.write("%s\n" % outScore)
outfile.close()



