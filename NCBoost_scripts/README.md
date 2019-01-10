# NCBoost scripts

Three main bash scripts are used for downloading and preparing the input data and scoring of variants.  
NOTE: The bash scripts have been designed to work on Linux systems.  

## 1- Downloading the features using get_selection_features.sh

The script calls wget and gsutil to download the different feature annotation files from the original sources in the NCBoost_features/ folder.

It may be run using the following command line from the root NCBoost folder:
```
./NCBoost_scripts/get_selection_features.sh TRUE 
```
The total size of the feature data files is about 450Go, so it might be long...
On non-linux systems, we recommend to manually download in the NCBoost_features/ folder the files indicated [here](https://github.com/RausellLab/NCBoost/blob/master/NCBoost_features/README.md).

## 2- Preparing the feature files using prepare_feature_databases.sh
To make the annotation more efficient, most of the feature mapping will be performed using tabix or pytabix. This script will index the files that are not indexed yet.  
It indexes the bed and vcf files containing CDTS, mean Het and mean DAF and the three scores from the 1000 Genome Project.

You may run it with the following command line, from the main NCBoost folder:
```
./NCBoost_scripts/prepare_feature_databases.sh
```

## 3- Scoring using NCBoost_annotate.sh
This last bash script will perform the annotation and scoring of the [variants](https://github.com/RausellLab/NCBoost#variant-input-format).

You may run it with the following command line, from the root NCBoost folder.
```
./NCBoost_scripts/ncboost_annotate.sh inF outF 
```
This bash script will automatically call several scripts (in perl, python and R) to perform the following tasks:

### A- Assessment of the closest gene and annotation of the genomic region type
First, the bash script will launch the gene-based annotation of ANNOVAR to obtain the genomic region and the name of the closest gene for each variant, using the default gene position database downloaded along ANNOVAR standard distribution (with a -splicing_threshold of 10bp).

`perl annotate_variation.pl -out outF -build hg19 -splicing_threshold 10 inF humandb/`

ANNOVAR gene-based annotation adds two columns: the first giving the type of genomic region overlapping the variant relative to the closest gene (e.g. intron, UTR'5, etc), and the second containing the name of the closest gene. Note that at this step ANNOVAR can provide several genomic regions and genes for a given variant. This kind of conflict is explained below.

### B- Dealing with closest gene and type of region conflicts
As NCBoost uses gene-based features and the genomic region type, it is important that each variant is properly associated to a single gene. A variant is uniquely associated to one closest gene and region type. To that aim, `ncboost_annotate.sh` calls `clean_annovar.py`, which performs the following tasks:

First, exonic, splicing and ncRNA variants are filtered out.  
Then, mappings are handled on a region-specific basis as follows:

- Intergenic variants can be mapped to either:
  * no gene, in this case NCBoost will not be able to score the corresponding variant(s).
  * one gene, which will be kept,
  * two genes, with their respective distance to the variant. In this case, the closest protein-coding gene will be kept 

- Upstream and Downstream variants can be:
  * upstream or downstream one gene, with which it will be associated.
  * upstream and/or downstream two genes. In this case it will be associated with the closest protein-coding gene (according to our gene database) and its  corresponding genomic region.


- UTRs variants can be associated with:
  * one gene,
  * two genes (with the corresponding relative distance to the gene), and will eventually be associated with the closest.
  
- Intronic variants can be associated with:
  * one gene,
  * two genes. In this case the variant is associated with the first protein-coding gene reported by ANNOVAR.

### C- One-hot encoding of the non-coding region
The type of non-coding region overlapping each variant is one-hot-encoded, and will be used by NCBoost as a feature.

### D- Adding gene-level features
Then, it adds gene-level features, through the gene symbol of the gene associated with each variant. The gene-based features used here are all contained in the file [NCBoost_geneDB.tsv](https://github.com/RausellLab/NCBoost/blob/master/NCBoost_data/NCBoost_geneDB.tsv), available in `NCBoost_data/`.

### E- Adding position-level and window-level features

`get_cadd_features.py` is called to add the 11 CADD-extracted features (PhCons x3, PhyloP x3, GerpN, GerpS, Bstatistic, GC and CpG) to the variant.

`get_1000GP_possel_scores.py` is called to add Tajima's D, Fu&Li's D* and F*.

`get_1000GP_daf_het.py` is called to add the mean derived allele frequency and mean heterozygosity over 1000bp, excluding the DAF or Het at the variant locus.

`get_GnomAD_MAFs.py` is called to compute the mean MAF over 1000bp (excluding the MAF of the variant), using the general MAF from the 1000 Genome Project, the general MAF from GnomAD and 7 population-specific MAFs.

`get_CDTS.py` adds to each variant the context-dependant tolerance score corresponding to the 10bp bin containing the position of the variant.

### F- Scoring with NCBoost
Finally, the R script `NCBoost_score.R` is called to apply the bundle of 10 independently trained XGBoost models to the variants as explained [here](https://github.com/RausellLab/NCBoost/tree/master/NCBoost_data#bundle-of-10-independently-trained-xgboost-models).

## Contact
Please address comments and questions about NCBoost to barthelemy.caron@institutimagine.org and antonio.rausell@inserm.fr
