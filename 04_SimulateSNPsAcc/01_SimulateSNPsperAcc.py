#!/usr/bin/python
import math
import random
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
import h5py
from pygwas.core import genotype

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs to be taken for analysis", type="int")
inOptions.add_option("-a", "--acc_to_check", dest="acc", help="Index of the accession needed to be checked", type="int")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-t", "--file_num_snps", dest="file_num_snps", help="Output from the CalculateSNPseachAcc.py script", type="string")
inOptions.add_option("-k", "--kinship_matrix", dest="kinFile", help="Kinship matrix file", type="string")
inOptions.add_option("-s", "--avgKin", dest="outKin", help="Output average kinship matrix", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")

(options, args) = inOptions.parse_args()

GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

AccToCheck = options.acc
numSNPs = options.numsnps

sampleSNPs = numpy.sort(random.sample(numpy.where(GenotypeData_acc.snps[:,AccToCheck] == 1)[0], numSNPs))
ScoreList = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")
for i in range(0, len(sampleSNPs), 1000):
  subsampleSNPs = GenotypeData.snps[tuple(sampleSNPs[i:i+1000]),]
  subsampleSNPs[subsampleSNPs < 0] = 0
  ScoreList = ScoreList + numpy.sum(subsampleSNPs, axis = 0)

# Calculate the number of SNPs in all the accessions
TotnumberSNPs = open(options.file_num_snps, 'r').read().split("\n")

FinalScore = numpy.zeros(len(GenotypeData.accessions))
FinalScore = [(float(ScoreList[i])*float(ScoreList[i]))/(numSNPs * float(TotnumberSNPs[i].split("\t")[1])) for i in range(0, len(GenotypeData.accessions))]
TopHitAccs = numpy.where(FinalScore > (max(FinalScore) - 2*numpy.std(FinalScore)))[0]
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  numsnp = TotnumberSNPs[i].split("\t")[1]
  outfile.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), int(numsnp), FinalScore[i]))

outfile.close()

#KinshipMat = h5py.File(options.kinFile)['kinship']
#outavgKin = open(options.outKin, 'a')
if len(TopHitAccs) > 1:
# Get the positions where the ambiguous acc actually differ
  outrefScore=open(options.refScore,'w')
  AccSNPs = GenotypeData_acc.snps[:, TopHitAccs]
  totAccSNPs = numpy.copy(AccSNPs)
  max= len(TopHitAccs)
  AccSNPs[AccSNPs == -1] = max + 100
  sumNumpy = numpy.sum(AccSNPs, axis = 1)
  RefPosInd = numpy.where((sumNumpy > 0) & (sumNumpy < max))[0]

  totAccSNPs[totAccSNPs < 0] = 0
  totalAccSNPs = numpy.sum(totAccSNPs, axis = 0)
  refScoreList = numpy.zeros(len(TopHitAccs))
  matSampleSNPs = RefPosInd[numpy.where(numpy.in1d(RefPosInd, sampleSNPs))]
  numMatched = len(matSampleSNPs)
  for i in range(0, len(matSampleSNPs), 1000):
    subsampleSNPs = totAccSNPs[tuple(matSampleSNPs[i:i+1000]),]
    subsampleSNPs[subsampleSNPs < 0] = 0
    refScoreList = refScoreList + numpy.sum(subsampleSNPs, axis = 0)

  for i in range(0, len(TopHitAccs)):
    outScore = (refScoreList[i]*refScoreList[i])/(numMatched*totalAccSNPs[i])   
    outrefScore.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], int(refScoreList[i]), totalAccSNPs[i], outScore))

  outrefScore.close()
#  subKin = KinshipMat[tuple(TopHitAccs),:][:,tuple(TopHitAccs)]
#  getKin = subKin[numpy.where(subKin != 1)]
#  KinScore = numpy.mean(getKin)
#  minKin = numpy.minimum(getKin)
#  maxKin = numpy.maximum(getKin)
#  outavgKin.write("%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[AccToCheck], len(TopHitAccs), minKin, maxKin, KinScore))

#outavgKin.close()

