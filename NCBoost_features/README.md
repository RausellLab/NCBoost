# Source feature files

This readme file reports the URLs of the files downloaded from various sources from which a number of purifying selection features were obtained.  
Files can be downloaded using the script NCBoost_scripts/get_selection_features.sh, from the NCBoost/ root folder.  
Please follow the download instructions and data processing described [here](https://github.com/RausellLab/NCBoost#2-Download-and-processing-of-the-feature-files-mined-by-NCBoost).

## Sources
### CADD v1.3 [[1]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references):
CADD is a variant classifier exploiting various conservation and epigenetic features.  
The annotated whole genome (as a tabix indexed file, ~350GO) can be downloaded here:  
`http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz`,  
with the corresponding tabix index file (~70Mo):  
`http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz.tbi`.  
From CADD v1.3 we extracted the PhyloP, PhastCons and Gerp scores, as well as the B statistic and both GC and CpG contents.

### Human population-specific natural selection scores
Human population-specific natural selection scores (Tajima's D, Fu & Li's D* and Fu & Li's F* percentiles for CEU, CHB and YRI populations) from the 1000 Genome Selection Browser [[2]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references) (9 files, ~20Mo each):
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisD_CEU.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisD_CHB.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisD_YRI.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisF_CEU.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisF_CHB.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisF_YRI.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/TajimasD_CEU.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/TajimasD_CHB.whole_genome.pvalues.gz`
`http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/TajimasD_YRI.whole_genome.pvalues.gz`
    
### Context-dependent tolerance scores (CDTS) [[3]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references) (~4.6Go):
`http://www.hli-opendata.com/noncoding/Pipeline/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed.gz`
    
    
### Mean Derived Allele Frequency and Mean Heterozygosity
Mean Derived Allele Frequency and Mean Heterozygosity of variants in a 1kb window region centered on the SNV and calculated from the 1000 Genomes Project (excluding the query variant) [[4]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references) (~450Mo each):
`ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/source_data/1kg/daf.bed.gz`
`ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/source_data/1kg/het_rates.bed.gz`


### GnomAD MAF [[5]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references)
GnomAD Minor Allele Frequency (~91Go):
`http://gnomad.broadinstitute.org/downloads` , VCF and coverage files, Genome Data (indexed vcf files and .tbi files for each chromosome).


## Feature details
* PhastCons [[6]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references) scores for three multi-species alignment (Vertebrates, Mammals and Primates, excluding human).
  
   PhastCons scores estimate the probability that the locus is contained in a conserved element.


* PhyloP [[7]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references) scores for three multi-species alignment (Vertebrates, Mammals and Primates, excluding human).
  
   PhyloP scores measure neutral evolution at individual sites. The score corresponds to the -log p-value of the null hypothesis of neutral evolution. 
   Positive values (up to 3) represent purifying selection, while negative values (up to -14) represent acceleration.


* GerpN and GerpS [[8]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references) single-nucleotide scores from mammalian alignments (excluding human).
  
   GerpN and GerpS single-nucleotide scores assess respectively the neutral substitution rate and the rejected substitution rate of the locus.
   A high GerpN value indicates high homology of the locus across species. 
   Positive values of GerpS indicate a deficit in substitutions, while negative values convey a substitution surplus.


* Tajima's D [[9]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references), Fu & Li's D* and Fu & Li's F* [[10]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references).
  
   Tajima's D is a neutrality test comparing estimates of the number of segregating sites and the mean pair-wise difference between sequences.
   Fu & Li's D* is a neutrality test comparing the number of singletons with the total number of nucleotide variants within a region.  
   Fu & Li's F* is a neutrality test comparing the number of singletons with the average number of nucleotide differences between pairs of sequences.  
   The three tests were performed within 3 populations of the 1000 Genome Project phase 1 data, producing population-specific scores: Yoruba in Ibadan, Nigeria (YRI), Han Chinese in Beijing, China (CHB) and Utah Residents with Northern and Western European Ancestry (CEU).  
   Negative logarithmic percentiles corresponding to each of these score (more precisions [there](http://hsb.upf.edu/?page_id=594)) were used with values ranging from 0 (indicating positive selection) to 6 (indicating purifying selection).  


* The background selection score - B statistic, [[11]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_features#references).
  
   It indicates the expected fraction of neutral diversity that is present at a site.  
   B statistic values close to 0 represent nearly complete removal of diversity as a result of selection and values near 1 indicate no conservation.  
   B-statistic is based on human single nucleotide polymorphism (SNP) data from Perlegen Sciences, HapMap phase II, the SeattleSNPs NHLBI Program for Genomic Applications and the NIEHS Environmental Genome Project. 


* Context-dependent tolerance scores (CDTS) for 10bp bins of the human genome computed on 15496 unrelated individuals from the gnomAD consortium. 
  
   The CDTS represents the difference between observed and expected variations in humans.  
   The expected variation is computed genome-wide for each nucleotide as the probability of variation of each nucleotide depending on its heptanucleotide context.  
   Low CDTS indicate loci intolerant to variation.


* Mean heterozygosity and mean derived allele frequency. 
  
   The Mean_Het and Mean_MAF of variants in a 1kb window region centered around the SNV and calculated from the 1000 Genomes Project (excluding the query variant).  
   Mean minor allele frequency (MAF) of variants in 1kb window were calculated from GnomAD genome data, excluding the query variant from the calculation.  
   Mean MAF was assessed for the global population and for population-specific frequencies: Africans and African Americans (AFR), Admixed Americans (AMR), East Asians (EAS), Finnish (FIN), Non-Finnish Europeans (NFE), Ashkenazi Jewish (ASJ) and Other populations (OTH).  
   Additionally, we extracted mean MAF of variants in a 1kb window calculated from the 1000 Genomes Project (excluding the query variant).


## References
1: Kircher et al. (2014). A general framework for estimating the relative pathogenicity of human genetic variants. Nat. Genet. 46, 310-315.

2: Pybus et al. (2014). 1000 Genomes Selection Browser 1.0: a genome browser dedicated to signatures of natural selection in modern humans. Nucleic Acids Res. 42, D903-D909.

3: di Iulio et al. (2018). The human noncoding genome defined by genetic diversity. Nat. Genet. 50, 333-337.

4: Ritchie et al. (2014). Functional annotation of noncoding sequence variants. Nat. Methods 11, 294-296.

5: Lek et al. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature 536, 285-291.

6: Hubisz et al. (2011). PHAST and RPHAST: phylogenetic analysis with space/time models. Brief. Bioinform. 12, 41-51.

7: Pollard et al. (2010). Detection of nonneutral substitution rates on mammalian phylogenies. Genome Res. 20, 110-121.

8: Davydov et al. (2010). Identifying a High Fraction of the Human Genome to be under Selective Constraint Using GERP++. PLoS Comput. Biol. 6, e1001025.

9: Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics 123, 585-595.

10: Fu and Li. (1993). Statistical tests of neutrality of mutations. Genetics 133, 693-709.

11: McVicker et al. (2009). Widespread Genomic Signatures of Natural Selection in Hominid Evolution. PLoS Genet. 5, e1000471.

## Contact
Please address comments and questions about NCBoost to barthelemy.caron@institutimagine.org and antonio.rausell@inserm.fr
