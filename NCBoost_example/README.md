# NCBoost example

The use of NCBoost is detailed here through a toy example.

## 1 - Getting an example input file

The [input file](https://github.com/RausellLab/NCBoost/blob/master/NCBoost_example/toy_variants.tsv) reports 100 common variants from the 1000 Genome Project in the format required by NCBoost.  


The head of the file looks like this:
```
chr	start		end		ref	alt	rsid
20	56092402	56092402	G	A	rs1985098
19	4302793		4302793		G	A	rs4807577
20	45035213	45035213	T	C	rs401379
12	57635257	57635257	C	T	rs7299890
12	49075984	49075984	T	G	rs61942053
13	41195992	41195992	T	C	rs4941987
```

Of note, the _rsid_ column is not mandatory.

## 2 - Running NCBoost

From the folder `NCBoost/`, run the script `ncboost_annotate.sh`, (available in `NCBoost_scripts/`), indicating as arguments the paths to the input and output files:
```
./NCBoost_scripts/ncboost_annotate.sh NCBoost_example/toy_variants.tsv NCBoost_example/toy_variants_annotated.tsv
```
The script will provide the feature annotation of variants and the NCBoost score
The terminal output should look as follows:

```
>./NCBoost_scripts/ncboost_annotate.sh NCBoost_example/toy_variants.tsv NCBoost_example/toy_variants_annotated.tsv


sorting data



Launching ANNOVAR

NOTICE: The --geneanno operation is set to ON by default
NOTICE: Output files were written to NCBoost_example/toy_variants.variant_function, NCBoost_example/toy_variants.exonic_variant_function
NOTICE: Reading gene annotation from annovar/humandb/hg19_refGene.txt ... Done with 52068 transcripts (including 11837 without coding sequence annotation) for 26464 unique genes
NOTICE: Processing next batch with 100 unique variants in 100 input lines
NOTICE: Variants with invalid input format were written to NCBoost_example/toy_variants.invalid_input


Cleanning ANNOVAR annotations, adding gene-based and context-awareness features

100%|#################################################################| 100/100 [00:00<00:00, 199.07it/s]
0 exonic/splicing/ncRNA associated positions have been removed
100 positions have been retained, among 100 total lines


Adding CADD features

100%|##################################################################| 101/101 [00:09<00:00, 10.98it/s]


Adding 1000 Genome Project positive selection scores

100%|#################################################################| 101/101 [00:00<00:00, 160.71it/s]
2 variants were not annotated due to lacking position in the database - scores are only available for autosomes


Adding 1000 Genome Project mean heterozygosity and DAF

100%|#################################################################| 101/101 [00:00<00:00, 438.39it/s]


Adding GnomAD MAFs

100%|##################################################################| 101/101 [00:11<00:00,  9.15it/s]


Adding CDTS score

100%|#################################################################| 101/101 [00:00<00:00, 277.14it/s]


Running NCBoost

	 NCBoost progressing

  |======================================================================| 100%

 50 var/sec 

```

## Contact
Please address comments and questions about NCBoost to barthelemy.caron@institutimagine.org and antonio.rausell@inserm.fr

