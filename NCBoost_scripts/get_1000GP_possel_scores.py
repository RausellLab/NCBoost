# coding: utf-8
# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 
# This script adds Fu&Li's D*, Fu&Li's F* and Tajimas'D for three populations.

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

tajimasD_YRI_data = tabix.open('NCBoost_features/TajimasD_YRI.whole_genome.pvalues.vcf.gz')
tajimasD_CEU_data = tabix.open('NCBoost_features/TajimasD_CEU.whole_genome.pvalues.vcf.gz')
tajimasD_CHB_data = tabix.open('NCBoost_features/TajimasD_CHB.whole_genome.pvalues.vcf.gz')

FuLisD_YRI_data = tabix.open('NCBoost_features/FuLisD_YRI.whole_genome.pvalues.vcf.gz')
FuLisD_CEU_data = tabix.open('NCBoost_features/FuLisD_CEU.whole_genome.pvalues.vcf.gz')
FuLisD_CHB_data = tabix.open('NCBoost_features/FuLisD_CHB.whole_genome.pvalues.vcf.gz')

FuLisF_YRI_data = tabix.open('NCBoost_features/FuLisF_YRI.whole_genome.pvalues.vcf.gz')
FuLisF_CEU_data = tabix.open('NCBoost_features/FuLisF_CEU.whole_genome.pvalues.vcf.gz')
FuLisF_CHB_data = tabix.open('NCBoost_features/FuLisF_CHB.whole_genome.pvalues.vcf.gz')

absent_index = 0

PosSel_annotations = ['TajimasD_YRI_pvalue', 'TajimasD_CEU_pvalue', 'TajimasD_CHB_pvalue', 'FuLisD_YRI_pvalue', 'FuLisD_CEU_pvalue', 'FuLisD_CHB_pvalue', 'FuLisF_YRI_pvalue', 'FuLisF_CEU_pvalue', 'FuLisF_CHB_pvalue']

with open(inF, 'r') as f, open(outF, 'w') as o:
    for index, l in enumerate(tqdm(f, total=get_num_lines(inF))):
        if l.startswith('chr'):
            old_header = l.rstrip().split('\t')
            standard_header = old_header[0:32]
            extra_header = old_header[32:len(old_header)]
            o.write("%s\n" % ('\t'.join(str(x) for x in (standard_header + PosSel_annotations + extra_header))))
            continue
        t = l.rstrip().split('\t')
        chrom = t[0]
        pos = int(t[1])
        standard = t[0:32]
        extra = t[32:len(t)]
        #get tajimasD YRI for the variant
        try:
            values = tajimasD_YRI_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
            absent_index += 1
            # print("index error handled, pursuing...", chrom, pos)
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        TajimasD_YRI_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'

    #get tajimasD CEU for the variant
        try:
            values = tajimasD_CEU_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        TajimasD_CEU_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'

    #get tajimasD CHB for the variant
        try:
            values = tajimasD_CHB_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        TajimasD_CHB_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'

    #get FuLisD YRI for the variant
        try:
            values = FuLisD_YRI_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        FuLisD_YRI_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'

    #get FuLisD CEU for the variant
        try:
            values = FuLisD_CEU_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        FuLisD_CEU_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'

    #get FuLisD CHB for the variant
        try:
            values = FuLisD_CHB_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        FuLisD_CHB_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'
    
    #get FuLisF YRI for the variant
        try:
            values = FuLisF_YRI_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        FuLisF_YRI_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'
   
    #get FuLisF CEU for the variant
        try:
            values = FuLisF_CEU_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        FuLisF_CEU_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'
        
    #get FuLisF CHB for the variant
        try:
            values = FuLisF_CHB_data.querys("%s:%s-%s" % (chrom, pos-3000, pos+3000))
        except tabix.TabixError:
            values = []
        data=[d for d in values]
        if len(data) != 0:
            start_pos = [int(d[1]) for d in data]
            reg_pos = [loc for loc in start_pos if loc <= pos]
            if len(reg_pos) > 1:
                reg_pos = reg_pos[1]
            if type(reg_pos) == int:
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            elif type(reg_pos) == list and len(reg_pos) > 0:
                reg_pos = reg_pos[0]
                scores = [d[3:5] for d in data if int(d[1]) == reg_pos][0]
            else:
                scores = []
        else:
            scores = []
        FuLisF_CHB_pvalue = float(scores[1]) if len(scores) == 2 else 'NA'
        PosSel_current = [TajimasD_YRI_pvalue, TajimasD_CEU_pvalue, TajimasD_CHB_pvalue, FuLisD_YRI_pvalue, FuLisD_CEU_pvalue, FuLisD_CHB_pvalue, FuLisF_YRI_pvalue, FuLisF_CEU_pvalue, FuLisF_CHB_pvalue]
        o.write("%s\n" % ('\t'.join(str(x) for x in (standard + PosSel_current + extra))))

print('{} variants were not annotated due to lacking position in the database - scores are only available for autosomes'.format(absent_index))
