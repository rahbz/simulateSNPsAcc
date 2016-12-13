#!/usr/bin/python
import math
import random
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import logging
import numpy
import numpy.ma
import h5py
from pygwas.core import genotype
import scipy
import statsmodels.stats.proportion as prop

#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs to be taken for analysis", type="int")
inOptions.add_option("-a", "--acc_to_check", dest="acc", help="Index of the accession needed to be checked", type="int")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-r", "--refScore", dest="refScore", help="Output for refined score", type="string")
inOptions.add_option("-t", "--pvalThres", dest="pVal", help="Threshold on P-value", default=0.95, type="float")
inOptions.add_option("-s", "--error_rate", dest="error", help="Maximum score which is considered to be for top hit accession", default=0.99, type="float")

(options, args) = inOptions.parse_args()
logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

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
  t1001SNPs = GenotypeData.snps[tuple(sampleSNPs[i:i+chunk_size]),]
  sam1kSNPs = t1001SNPs[:, AccToCheck]
  samSNPs = numpy.reshape(numpy.repeat(sam1kSNPs, num_lines), (len(sam1kSNPs), num_lines))
  ScoreList = ScoreList + numpy.sum(t1001SNPs == samSNPs, axis=0)
  NumInfoSites = NumInfoSites + len(sampleSNPs[i:i+chunk_size]) - numpy.sum(numpy.ma.masked_less(t1001SNPs, 0).mask, axis = 0)
  logging.info("Done analysing %s SNPs", i+chunk_size)


pVals = [prop.proportions_ztest(ScoreList[i], NumInfoSites[i], value=options.error, alternative="smaller")[1] for i in range(num_lines)]
pVals = numpy.array(pVals).astype("float")
outfile = open(options.outFile, 'w')
for i in range(0, len(GenotypeData.accessions)):
  outfile.write("%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[i], int(ScoreList[i]), NumInfoSites[i], float(ScoreList[i]) / NumInfoSites[i], pVals[i]))

outfile.close()
logging.info("Written scores into the output file")

TopHitAccs = numpy.where(pVals > options.pVal)[0]
if len(TopHitAccs) > 1:
## Get the positions where the ambiguous acc actually differ
  outrefScore=open(options.refScore,'w')
  AccSNPs = GenotypeData_acc.snps[:, TopHitAccs]
  AccReq = numpy.where(TopHitAccs == options.acc)[0]
  max = len(TopHitAccs)
  AccSNPs[AccSNPs == -1] = max + 3
  sumNumpy = numpy.sum(AccSNPs, axis = 1)
  RefPosInd = numpy.where((sumNumpy > 0) & (sumNumpy < max))[0]
  infoSites = RefPosInd[numpy.in1d(RefPosInd, sampleSNPs)]
  AccGTs = AccSNPs[:, AccReq]
  samGTs = AccGTs[infoSites]
  samSNPs = numpy.reshape(numpy.repeat(samGTs, max), (len(infoSites), max))
  refScoreList = numpy.sum(AccSNPs[infoSites, :]  == samSNPs, axis=0)
  for i in range(0, max):
    rpval = prop.proportions_ztest(refScoreList[i], len(infoSites), value=options.error, alternative="smaller")[1]
    outrefScore.write("%s\t%s\t%s\t%s\t%s\n" % (GenotypeData.accessions[TopHitAccs[i]], refScoreList[i], len(infoSites), float(refScoreList[i])/len(infoSites), rpval))
  outrefScore.close()
