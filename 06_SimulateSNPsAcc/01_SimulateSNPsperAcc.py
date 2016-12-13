#!/usr/bin/python
import math
import random
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import numpy
import numpy.ma
import h5py
from pygwas.core import genotype

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs to be taken for analysis", type="int")
inOptions.add_option("-a", "--acc_to_check", dest="acc", help="Index of the accession needed to be checked", type="int")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-k", "--kinship_matrix", dest="kinFile", help="Kinship matrix file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")

(options, args) = inOptions.parse_args()

GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

AccToCheck = options.acc
numSNPs = options.numsnps
chunk_size = 1000

sampleSNPs = numpy.sort(random.sample(numpy.where(GenotypeData_acc.snps[:,AccToCheck] >= 0)[0], numSNPs))
ScoreList = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")
NumInfoSites = numpy.zeros(len(GenotypeData.accessions), dtype="uint32")

num_lines = len(GenotypeData.accessions)
for i in range(0, len(sampleSNPs), chunk_size):
  sam1kSNPs = GenotypeData_acc.snps[:,AccToCheck][sampleSNPs[i:i+chunk_size]]
  samSNPs = numpy.reshape(numpy.repeat(sam1kSNPs, num_lines), (len(sam1kSNPs), num_lines))
  t1001SNPs = GenotypeData.snps[tuple(sampleSNPs[i:i+chunk_size]),]
  tempBool = t1001SNPs == samSNPs
  tempBool = tempBool.astype(int)
  tempScore  = numpy.array(numpy.sum(tempBool, axis=0), dtype="uint16")
  NumInfoSites = NumInfoSites + len(sampleSNPs[i:i+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask, axis = 0)
  ScoreList = ScoreList + tempScore

FinalScore = numpy.zeros(len(GenotypeData.accessions))
FinalScore = [float(ScoreList[i]) / float(NumInfoSites[i]) for i in range(0, len(GenotypeData.accessions))]

TopHitAccs = numpy.where(FinalScore > (max(FinalScore) - 2*numpy.std(FinalScore)))[0]
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  outfile.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], FinalScore[i]))

outfile.close()

#KinshipMat = h5py.File(options.kinFile)['kinship']
#outavgKin = open(options.outKin, 'a')
#if len(TopHitAccs) > 1:
## Get the positions where the ambiguous acc actually differ
#  outrefScore=open(options.refScore,'w')
#  AccSNPs = GenotypeData_acc.snps[:, TopHitAccs]
#  totAccSNPs = numpy.copy(AccSNPs)
#  max= len(TopHitAccs)
#  AccSNPs[AccSNPs == -1] = max + 100
#  sumNumpy = numpy.sum(AccSNPs, axis = 1)
#  RefPosInd = numpy.where((sumNumpy > 0) & (sumNumpy < max))[0]
#
#  totAccSNPs[totAccSNPs < 0] = 0
#  totalAccSNPs = numpy.sum(totAccSNPs, axis = 0)
#  refScoreList = numpy.zeros(len(TopHitAccs))
#  matSampleSNPs = RefPosInd[numpy.where(numpy.in1d(RefPosInd, sampleSNPs))]
#  numMatched = len(matSampleSNPs)
#  for i in range(0, len(matSampleSNPs), 1000):
#    subsampleSNPs = totAccSNPs[tuple(matSampleSNPs[i:i+1000]),]
#    subsampleSNPs[subsampleSNPs < 0] = 0
#    refScoreList = refScoreList + numpy.sum(subsampleSNPs, axis = 0)
#
#  for i in range(0, len(TopHitAccs)):
#    outScore = (refScoreList[i]*refScoreList[i])/(numMatched*totalAccSNPs[i])   
#    outrefScore.write("%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], int(refScoreList[i]), totalAccSNPs[i], outScore))
#
#  outrefScore.close()
#  subKin = KinshipMat[tuple(TopHitAccs),:][:,tuple(TopHitAccs)]
#  getKin = subKin[numpy.where(subKin != 1)]
#  KinScore = numpy.mean(getKin)
#  minKin = numpy.minimum(getKin)
#  maxKin = numpy.maximum(getKin)
#  outavgKin.write("%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[AccToCheck], len(TopHitAccs), minKin, maxKin, KinScore))

#outavgKin.close()

