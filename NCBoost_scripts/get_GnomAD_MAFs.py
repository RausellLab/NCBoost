# coding: utf-8
# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 
# This script compute and adds the 8 mean MAFs from GnomAD and the mean MAF from the 1000GP.

import sys, vcf, tabix, numpy, mmap
from decimal import Decimal
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
    

def getPopulationL(ac_pop, an_pop):
    acL = list(map(float, ac_pop))
    # if an_pop is int:
    # 	an_pop = [an_pop]
    # an = numpy.sum(map(float, an_pop))
    an = float(an_pop)
    if acL != [] and acL != [None]:
        rac = an - numpy.sum(acL)
        acL.append(rac)
        mac = sorted(acL, reverse=True)[1]
        af = mac / an if an != 0 else None
    return af 

def listFromRegion(chrom, start, end, gnomadAF):
    header = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|ExAC_MAF|ExAC_Adj_MAF|ExAC_AFR_MAF|ExAC_AMR_MAF|ExAC_EAS_MAF|ExAC_FIN_MAF|ExAC_NFE_MAF|ExAC_OTH_MAF|ExAC_SAS_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info".split('|')
    gmaf_idx = header.index('GMAF')
    gmafL = []
    mafL = []
    maf_afrL = []
    maf_amrL = []
    maf_asjL = []
    maf_easL = []
    maf_finL = []
    maf_nfeL = []
    maf_othL = []	
    try:
        region = gnomadAF.fetch(chrom, start, end)  # pyvcf fetch is 0 based
    except ValueError:
        print("region position causing problems, variant not treated, pursuing...", chrom, start, end)
        region = []
    for r in region:
        afl = r.INFO['AF']
        afl = list(filter(lambda a:a != None, afl))
        if afl != []:
            raf = 1 - numpy.sum(afl)
            afl.append(raf)
            maf = sorted(afl, reverse=True)[1]
            mafL.append(maf)
            annotation = [None if a =='' else a for a in r.INFO['CSQ'][0].split('|')]
            gmaf = float(annotation[gmaf_idx].split('&')[0].split(':')[1]) if annotation[gmaf_idx] != None else None
            if gmaf != None: gmafL.append(gmaf)
        af_afr = getPopulationL(r.INFO['AC_AFR'], r.INFO['AN_AFR'])
        if af_afr != None: maf_afrL.append(af_afr)
        af_amr = getPopulationL(r.INFO['AC_AMR'], r.INFO['AN_AMR'])
        if af_amr != None: maf_amrL.append(af_amr)
        af_asj = getPopulationL(r.INFO['AC_ASJ'], r.INFO['AN_ASJ'])
        if af_asj != None: maf_asjL.append(af_asj)
        af_eas = getPopulationL(r.INFO['AC_EAS'], r.INFO['AN_EAS'])
        if af_eas != None: maf_easL.append(af_eas)
        af_fin = getPopulationL(r.INFO['AC_FIN'], r.INFO['AN_FIN'])
        if af_fin != None: maf_finL.append(af_fin)
        af_nfe = getPopulationL(r.INFO['AC_NFE'], r.INFO['AN_NFE'])
        if af_nfe != None: maf_nfeL.append(af_nfe)
        af_oth = getPopulationL(r.INFO['AC_OTH'], r.INFO['AN_OTH'])
        if af_oth != None: maf_othL.append(af_oth)
    return gmafL, mafL, maf_afrL, maf_amrL, maf_asjL, maf_easL, maf_finL, maf_nfeL, maf_othL

gnomad_annotations = ['meanMAF1000G', 'meanMAFGnomAD', 'meanMAF_AFRGnomAD', 'meanMAF_AMRGnomAD', 'meanMAF_ASJGnomAD', 'meanMAF_EASGnomAD', 'meanMAF_FINGnomAD', 'meanMAF_NFEGnomAD', 'meanMAF_OTHGnomAD']

with open(inF, 'r') as f, open(outF, 'w') as o:
    for index, l in enumerate(tqdm(f, total=get_num_lines(inF))): 
        if l.startswith('chr'):
            old_header = l.rstrip().split('\t')
            standard_header = old_header[0:43]
            extra_header = old_header[43:len(old_header)]
            o.write("%s\n" % ('\t'.join(str(x) for x in (standard_header + gnomad_annotations + extra_header))))
            continue
        t = l.rstrip().split('\t')
        chrom = t[0]
        pos = int(t[1])
        standard = t[0:43]
        extra = t[43:len(t)]
        if chrom != "Y": 
            gnomadAF = vcf.Reader(filename='NCBoost_features/gnomad.genomes.r2.0.2.sites.chr{}.vcf.bgz'.format(chrom), compressed=True)
            # af of the 1KB region
            # 500nt left part
            left_gmafL, left_mafL, left_maf_afrL, left_maf_amrL, left_maf_asjL, left_maf_easL, left_maf_finL, left_maf_nfeL, left_maf_othL = listFromRegion(chrom, pos-501, pos-1, gnomadAF)
            # 500nt right part
            right_gmafL, right_mafL, right_maf_afrL, right_maf_amrL, right_maf_asjL, right_maf_easL, right_maf_finL, right_maf_nfeL, right_maf_othL = listFromRegion(chrom, pos, pos+500, gnomadAF)
            gmafL = left_gmafL + right_gmafL
            mafL = left_mafL + right_mafL
            maf_afrL = left_maf_afrL + right_maf_afrL
            maf_amrL = left_maf_amrL + right_maf_amrL
            maf_asjL = left_maf_asjL + right_maf_asjL
            maf_easL = left_maf_easL + right_maf_easL
            maf_finL = left_maf_finL + right_maf_finL
            maf_nfeL = left_maf_nfeL + right_maf_nfeL
            maf_othL = left_maf_othL + right_maf_othL
            meanGMAF = '%.2E' % Decimal(numpy.mean(gmafL)) if gmafL != [] else 'NA'
            meanMAF = '%.2E' % Decimal(numpy.mean(mafL)) if mafL != [] else 'NA'
            meanMAF_AFR = '%.2E' % Decimal(numpy.mean(maf_afrL)) if maf_afrL != [] else 'NA'
            meanMAF_AMR = '%.2E' % Decimal(numpy.mean(maf_amrL)) if maf_amrL != [] else 'NA'
            meanMAF_ASJ = '%.2E' % Decimal(numpy.mean(maf_asjL)) if maf_asjL != [] else 'NA'
            meanMAF_EAS = '%.2E' % Decimal(numpy.mean(maf_easL)) if maf_easL != [] else 'NA'
            meanMAF_FIN = '%.2E' % Decimal(numpy.mean(maf_finL)) if maf_finL != [] else 'NA'
            meanMAF_NFE = '%.2E' % Decimal(numpy.mean(maf_nfeL)) if maf_nfeL != [] else 'NA'
            meanMAF_OTH = '%.2E' % Decimal(numpy.mean(maf_othL)) if maf_othL != [] else 'NA'
        else:
            meanGMAF = meanMAF = meanMAF_AFR = meanMAF_AMR = meanMAF_ASJ = meanMAF_EAS = meanMAF_FIN = meanMAF_NFE = meanMAF_OTH = 'NA'
        gnomad_current = [meanGMAF, meanMAF, meanMAF_AFR, meanMAF_AMR, meanMAF_ASJ, meanMAF_EAS, meanMAF_FIN, meanMAF_NFE, meanMAF_OTH]
        o.write("%s\n" % ('\t'.join(str(x) for x in (standard + gnomad_current + extra))))
