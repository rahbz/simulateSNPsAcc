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
inOptions.add_option("-s", "--score_acc", dest="scoreAcc", help="Score Acc text file", type="string")
inOptions.add_option("-k", "--kinship_file", dest="kinFile", help="HDF5 file for kinship matrix", type="string")

inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the average kinship score", type="string")


(options, args) = inOptions.parse_args()

KinshipMat = h5py.File(options.kinFile)['kinship']
KinAcc = h5py.File(options.kinFile)['accessions']

AccToCheck = options.acc
numSNPs = options.numsnps

# Calculate the number of SNPs in all the accessions
ScoreFile = open(options.scoreAcc, 'r').read().split("\n")
Scores = [float(ScoreFile[i].split("\t")[3]) for i in range(0, len(ScoreFile)) if len(ScoreFile[i].split("\t")) == 4]

AmbAcc = numpy.where(Scores > numpy.max(Scores) - (2*numpy.std(Scores)))[0]

outAvg = open(options.outFile, 'a')
if len(AmbAcc) > 1:
  subKin = KinshipMat[tuple(AmbAcc),:][:,tuple(AmbAcc)] 
  non1subKin = numpy.where(subKin != 1)
  minKin = numpy.min(subKin[non1subKin])
  maxKin = numpy.max(subKin[non1subKin])
  KinScore = numpy.mean(subKin[non1subKin])
  outAvg.write("%s\t%s\t%s\t%s\t%s\n" % (KinAcc[AccToCheck], len(AmbAcc), minKin, maxKin, KinScore))

outAvg.close()


