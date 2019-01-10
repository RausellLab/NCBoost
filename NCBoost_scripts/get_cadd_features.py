# coding: utf-8
# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 
# This script uses tabix to get the bstatistic, the interspecies conservation features and both GC and CpG from the CADD whole genome file.

import sys, tabix, mmap
from tqdm import tqdm
import numpy as np

argvL = sys.argv
inF = argvL[1]
outF = argvL[2]
path_to_cadd = argvL[3]

def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines

CADD_data = tabix.open('{}/whole_genome_SNVs_inclAnno.tsv.gz'.format(path_to_cadd))

cadd_annotations = ["GC", "CpG", "priPhCons", "mamPhCons", "verPhCons", "priPhyloP", "mamPhyloP", "verPhyloP", "GerpN", "GerpS", "bStatistic"]

feature_indexes = [13,14,18,19,20,21,22,23,24,25,28]



with open(inF, 'r') as f, open(outF, 'w') as fo:
	for index, l in enumerate(tqdm(f, total=get_num_lines(inF))):
		cadd_values = []
		if l.startswith('chr'):
		  old_header = l.rstrip().split('\t')
		  standard_header = old_header[0:21]
		  extra_header = old_header[21:len(old_header)]
		  fo.write("%s\n" % ('\t'.join(str(x) for x in (standard_header + cadd_annotations + extra_header))))
		  continue
		t = l.rstrip().split('\t')
		chrom = t[0]
		pos = int(t[1])
		ref = t[2]
		alt = t[3]
		standard = t[0:21]
		extra = t[21:len(t)]
		#get CADD features
		try:
			values = CADD_data.querys("%s:%s-%s" % (chrom, pos, pos))
		except tabix.TabixError:
		  values = []
		  print("index error handled, pursuing...", "SNV {}".format(index), chrom, pos)
		data=[d for d in values]
		if len(data) != 0:
		  data2=[l_in for l_in in data if alt == l_in[4]]
		  cadd_values=np.take(data2[0], [13,14,18,19,20,21,22,23,24,25,28])
		CADD = cadd_values if len(cadd_values) != 0 else np.array(['NA' for x in feature_indexes]).astype(str)
		fo.write("%s\n" % ('\t'.join(str(x) for x in (standard + list(CADD) + extra))))
