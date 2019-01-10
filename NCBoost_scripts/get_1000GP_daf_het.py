# coding: utf-8
# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 
# This script adds the mean Derived Allele Frequency and Heterozygosity over the 1kb surronding the variant,
# excluding the DAF/Het of the variant itself.

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
    
daf = tabix.open('NCBoost_features/daf.bed.gz')
het = tabix.open('NCBoost_features/het_rates.bed.gz')

daf_het_annotations = ['meanDaf1000G', 'meanHet1000G']
  
with open(inF, 'r') as f, open(outF, 'w') as o:
	for l in tqdm(f, total=get_num_lines(inF)): 
		if l.startswith('chr'):
			old_header = l.rstrip().split('\t')
			standard_header = old_header[0:41]
			extra_header = old_header[41:len(old_header)]
			o.write("%s\n" % ('\t'.join(str(x) for x in (standard_header + daf_het_annotations + extra_header))))
			continue
		t = l.rstrip().split('\t')
		chrom = t[0]
		pos = int(t[1])
		standard = t[0:41]
		extra = t[41:len(t)]
		# daf of the postion
		try:
		  c_d = daf.querys("chr%s:%s-%s" % (chrom, pos, pos))
		except IndexError:
		  c_d = []
		  print("index error handled, pursuing...", chrom, pos)
		except tabix.TabixError:
		  c_d = []
		  print("tabix error handled, pursuing...", chrom, pos)
		c_daf_L = sum([c for c in c_d], [])
		c_daf = float(c_daf_L[4]) if len(c_daf_L) == 5 else None
		# daf of the region
		try:
		  dafL = daf.querys("chr%s:%s-%s" % (chrom, pos-500, pos+500))
		except (IndexError, tabix.TabixError) as error:
		  dafL = []		  
		DafL = [float(d[4]) for d in dafL]
		# remove the position daf if is present in the region
		if c_daf != None: 
			DafL.remove(c_daf)
		meandDaf = numpy.mean(DafL) if DafL != [] else 'NA'
		# het of the postion
		try:
		  c_h = het.querys("chr%s:%s-%s" % (chrom, pos, pos))
		except (IndexError, tabix.TabixError) as error:
		  c_h = []
		c_h_L = sum([c for c in c_h], [])
		c_het = float(c_h_L[4]) if len(c_h_L) == 5 else None
		# het of the region
		try:
		  hetL = het.querys("chr%s:%s-%s" % (chrom, pos-500, pos+500))
		except (IndexError, tabix.TabixError) as error:
		  hetL = []
		HetL = [float(h[4]) for h in hetL]
		# remove the position het if is present in the region
		if c_het != None:
			HetL.remove(c_het)
		meanHet = numpy.mean(HetL) if HetL != [] else 'NA'
		daf_het_current = [meandDaf, meanHet]
		o.write("%s\n" % ('\t'.join(str(x) for x in (standard + daf_het_current + extra))))
