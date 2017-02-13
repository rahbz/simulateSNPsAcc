#!/usr/bin/python
import math
import random
from optparse import OptionParser
#These are the modules that are needed for this script
# module load numpy
# module use /net/gmi.oeaw.ac.at/software/shared/nordborg_common/modulefiles/
# module load pygwas

import logging
import numpy as np
import numpy.ma
import h5py
from pygwas.core import genotype
import scipy
from snpmatch.core import snpmatch


#__________________________________________
inOptions = OptionParser()
inOptions.add_option("-n", "--num_snps_vcf", dest="numsnps", help="Number of SNPs to be taken for analysis", type="int")
inOptions.add_option("-a", "--acc_to_check", dest="accID", help="Index of the accession needed to be checked", type="string")
inOptions.add_option("-d", "--hdf5_file", dest="hdf5File", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-e", "--hdf5_acc_file", dest="hdf5accFile", help="Path to SNP matrix given in binary hdf5 file", type="string")
inOptions.add_option("-o", "--output", dest="outFile", help="Output file with the probability scores", type="string")
inOptions.add_option("-s", "--error_rate", dest="error", help="Maximum score which is considered to be for top hit accession", default=0.001, type="float")

(options, args) = inOptions.parse_args()
logging.basicConfig(format='%(levelname)s:%(asctime)s:  %(message)s', level=logging.DEBUG)

GenotypeData_acc = genotype.load_hdf5_genotype_data(options.hdf5accFile)
GenotypeData = genotype.load_hdf5_genotype_data(options.hdf5File)

AccID = options.accID
try:
  AccToCheck = np.where(GenotypeData.accessions == AccID)[0][0]
except:
  snpmatch.die("accession is not present in the matrix!")
numSNPs = options.numsnps
chunk_size = 1000

logging.info("randomly choosing %s SNPs from accession %s", numSNPs, AccID)
sampleSNPs = np.sort(random.sample(np.where(GenotypeData_acc.snps[:,AccToCheck] >= 0)[0], numSNPs))
ScoreList = np.zeros(len(GenotypeData.accessions), dtype="uint32")
NumInfoSites = np.zeros(len(GenotypeData.accessions), dtype="uint32")
NumMatSNPs = 0
num_lines = len(GenotypeData.accessions)

for i in range(0, len(sampleSNPs), chunk_size):
  t1001SNPs = GenotypeData.snps[tuple(sampleSNPs[i:i+chunk_size]),]
  sam1kSNPs = t1001SNPs[:, AccToCheck]
  samSNPs = np.reshape(np.repeat(sam1kSNPs, num_lines), (len(sam1kSNPs), num_lines))
  ScoreList = ScoreList + np.sum(t1001SNPs == samSNPs, axis=0)
  NumInfoSites = NumInfoSites + len(sampleSNPs[i:i+chunk_size]) - np.sum(numpy.ma.masked_less(t1001SNPs, 0).mask, axis = 0)
  NumMatSNPs = NumMatSNPs + len(sam1kSNPs)
  if i/chunk_size % 50 == 0:
    logging.info("Done analysing %s SNPs", NumMatSNPs)

logging.info("writing data!")

ScoreList = np.array([random.uniform(ScoreList[i] - (ScoreList[i] * options.error), ScoreList[i]) for i in range(num_lines)], dtype=int)

snpmatch.print_out_table(options.outFile, GenotypeData.accessions, ScoreList, NumInfoSites, numSNPs, "NA")

logging.info("finished!")
