# NCBoost gene-based features and bundle of 10 XGBoost models

## File NCBoost_geneDB.tsv
The file [NCBoost_geneDB.tsv](https://github.com/RausellLab/NCBoost/blob/master/NCBoost_data/NCBoost_geneDB.tsv) contains the 18744 protein-coding genes for which at least one gene-based feature was retrieved. During the processing of the source data files, gene names were mapped to approved HGNC gene symbols. Missing values were imputed through the median value computed over all protein-coding genes. The following gene-based features are reported: 

### - Primate dn/ds ratios (slr_dnds)
Primate dn/ds ratios (i.e. the ratio between the number of nonsynonymous substitutions and the number of synonymous substitutions), taken from [[1]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references). Low dn/ds values reflect purifying selection, while high dn/ds values are indicative of positive selection.

### - The gene probability of loss-of-function intolerance (pLI)
The pLI estimates the depletion of rare and de novo protein-truncating variants as compared to the expectations drawn from a neutral model of de novo variation on ExAC exomes data. pLI values close to 1 reflect gene intolerant to heterozygous and homozygous loss-of-function mutations. pLI was extracted from [[2]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).

### - Gene Damage Index (GDI)
GDI is a gene-level metric of the mutational damage that has accumulated in the general population, based on CADD scores. High GDI values reflect highly damaged genes.  
GDI was extracted from [[3]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).

### - The Residual Variation Intolerance Score (RVIS_percentile)
RVIS percentile estimates the gene departure from the average number of common functional mutations in genes with a similar amount of mutational burden in human. High RVIS percentiles reflect genes highly tolerant to variation. RVIS_percentile was extracted from [[4]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).

### - The non-coding version of the RVIS score (ncRVIS)
ncRVIS measures the departure from the genome-wide average of the number of common variants found in the noncoding sequence of genes with a similar amount of noncoding mutational burden in human. Negative values of ncRVIS indicate a conserved proximal non-coding region, while positive values indicate a higher burden of SNVs than expected under neutrality. ncRVIS was extracted from [[5]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).

### - The average non-coding GERP (ncGERP)
ncGERP is the average GERP score [[6]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references) across the non-coding sequence of a gene. ncGERP was extracted from [[5]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).

Both in the case of ncRVIS and ncGERP, the non-coding sequence was defined in the original publication as the collection of 5'-UTR, 3'-UTR and an additional non-exonic 250 bp upstream of transcription start site (TSS).  

### - The gene age (gene_age)
The gene age estimates the gene time of origin based on the presence/absence of orthologs in the vertebrate phylogeny. It varies from 0 (oldest) to 12 (youngest, corresponding to human specific genes). It was obtained from [[7]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).

### - The number of human paralogs (familyMemberCount)
The number of human paralogs for each gene was collected from the OGEE database [[8]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references).


## gene_to_sample
This file contains the protein-coding gene list, reporting the associated chromosome, the genomic partition to which a gene was assigned (sample_index; see below), and a pathogenic_flag, indicating whether the gene was associated with a variant present in our curated set of pathogenic non-coding variants.


## Bundle of 10 independently trained XGBoost models
NCBoost training was performed with XGBoost, a machine learning technique based on gradient tree boosting (also known as gradient boosted regression tree). To train NCBoost, we first randomly split the complete list of protein-coding genes in 10 genomic partitions of equal size, with the same distribution across all chromosomes and keeping in each of them the same proportion of monogenic Mendelian disease genes presenting high-confidence pathogenic non-coding variants. Each pathogenic SNV was associated with a unique set of 10 "negative" variants, randomly sampled from a set of common human variants without clinical assertion and associated with genes within the same genomic partition. Random sampling of common variants was matched to the positive set to keep the same proportion of variants per type of region: intronic, 5'UTR, 3'UTR, upstream, downstream and intergenic regions. A maximum of one positive and one negative variant associated with the same gene was allowed, although no minimum per gene was required. For the training step, a maximum of one pathogenic non-coding variant was randomly sampled per gene. We then trained NCBoost as a bundle of 10 independently trained XGBoost models, consecutively excluding in each of them 1 of the 10 genomic partitions described above. We note here that each SNV received one single score -and not 10 different scores-, which was provided by the model that did not contain the genomic partition wherein the SNV overlapped. Further details are provided in the [NCBoost paper](https://rdcu.be/bmlxX).

The _.Rdata_ files provided here contain the 10 trained xgboost [[9]](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#references) models that collectively constitute the NCBoost bundle, excluding in each of them the corresponding Xth genomic partition.

## References
1: Rausell et al. (2014). Analysis of Stop-Gain and Frameshift Variants in Human Innate Immunity Genes. PLoS Comput Biol 10, e1003757.

2: Lek et al. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature 536, 285-291.

3: Itan et al. (2015). The human gene damage index as a gene-level approach to prioritizing exome variants. Proc. Natl. Acad. Sci. 112, 13615-13620.

4: Petrovski et al. (2013). Genic Intolerance to Functional Variation and the Interpretation of Personal Genomes. PLoS Genet 9, e1003709.

5: Petrovski et al. (2015). The Intolerance of Regulatory Sequence to Genetic Variation Predicts Gene Dosage Sensitivity. PLOS Genet. 11, e1005492.

6: Davydov et al. (2010). Identifying a High Fraction of the Human Genome to be under Selective Constraint Using GERP++. PLoS Comput. Biol. 6, e1001025.

7: Popadin et al. (2014). Gene Age Predicts the Strength of Purifying Selection Acting on Gene Expression Variation in Humans. Am. J. Hum. Genet. 95, 660-674.

8: Chen et al. (2017). OGEE v2: an update of the online gene essentiality database with special focus on differentially essential genes in human cancer cell lines. Nucleic Acids Res. 45, D940-D944.

9: Chen and Guestrin. (2016). XGBoost: A Scalable Tree Boosting System. (ACM Press), pp. 785-794.

## Contact
Please address comments and questions about NCBoost to barthelemy.caron@institutimagine.org and antonio.rausell@inserm.fr


