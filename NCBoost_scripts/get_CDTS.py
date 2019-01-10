# coding: utf-8
# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 
# This script adds the CDTS score of the 10bp bin containing the variant.

import sys, tabix, numpy, mmap
from tqdm import tqdm

argvL = sys.argv
inF = argvL[1]
outF = argvL[2]

def get_num_lines(file_path):
    fp = open(file_path, "r+")
    buf = mmap.mmap(fp.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


CDTS_data = tabix.open('NCBoost_features/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed.gz')

cdts_annotations = ['CDTS']

with open(inF, 'r') as f, open(outF, 'w') as o:
	for index, l in enumerate(tqdm(f, total=get_num_lines(inF))):
		if l.startswith('chr'):
			old_header = l.rstrip().split('\t')
			standard_header = old_header[0:52]
			extra_header = old_header[52:len(old_header)]
			o.write("%s\n" % ('\t'.join(str(x) for x in (standard_header + cdts_annotations + extra_header))))
			continue
		t = l.rstrip().split('\t')
		chrom = t[0]
		pos = int(t[1])
		standard = t[0:52]
		extra = t[52:len(t)]
		#get CDTS for the variant's corresponding 10bp bin
		try:
			values = CDTS_data.querys("%s:%s-%s" % (chrom, pos-6, pos+6))
		except tabix.TabixError:
			values = []
			print("index error handled, pursuing...", "SNV n {}, chr {} pos {}".format(index, chrom, pos))
		data=[d for d in values]
		if len(data) != 0:
			start_pos = [int(d[1]) for d in data]
			reg_pos = [loc for loc in start_pos if loc <= pos]
			if len(reg_pos) > 1:
			  reg_pos = reg_pos[1]
			if type(reg_pos) == int:
			  scores = [d[5] for d in data if int(d[1]) == reg_pos][0]
			elif type(reg_pos) == list and len(reg_pos) > 0:
			  reg_pos = reg_pos[0]
			  scores = [d[5] for d in data if int(d[1]) == reg_pos][0]
			else:
			  scores = []
		else:
			scores = []
		CDTS = float(scores) if len(scores) != 0 else 'NA'
		CDTS_current = [CDTS]
		o.write("%s\n" % ('\t'.join(str(x) for x in (standard + CDTS_current + extra))))
