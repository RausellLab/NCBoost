# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 

folder='NCBoost_features'

### GWAVA - mean DAF, bed file, 0-based
gzip -d $folder/daf.bed.gz
bgzip $folder/daf.bed
tabix -p bed $folder/daf.bed.gz

### GWAVA - mean HET, bed file, 0-based
gzip -d $folder/het_rates.bed.gz
bgzip $folder/het_rates.bed
tabix -p bed $folder/het_rates.bed.gz


### 1000GP FuLi's D* and F*, vcf, 1-based
## FuLisD
for data in CEU CHB YRI; do
  gzip -d $folder/FuLisD_$data.whole_genome.pvalues.gz
  cut -d" " -f2- $folder/FuLisD_$data.whole_genome.pvalues > $folder/FuLisD_$data.whole_genome.pvalues_tmp.vcf
  awk -v OFS="\t" '$1=$1' $folder/FuLisD_$data.whole_genome.pvalues_tmp.vcf > $folder/FuLisD_$data.whole_genome.pvalues.vcf
  bgzip $folder/FuLisD_$data.whole_genome.pvalues.vcf
  tabix -p vcf -S 1 -f $folder/FuLisD_$data.whole_genome.pvalues.vcf.gz
done

## FuLisF
for data in CEU CHB YRI; do
  gzip -d $folder/FuLisF_$data.whole_genome.pvalues.gz
  cut -d" " -f2- $folder/FuLisF_$data.whole_genome.pvalues > $folder/FuLisF_$data.whole_genome.pvalues_tmp.vcf
  awk -v OFS="\t" '$1=$1' $folder/FuLisF_$data.whole_genome.pvalues_tmp.vcf > $folder/FuLisF_$data.whole_genome.pvalues.vcf
  bgzip $folder/FuLisF_$data.whole_genome.pvalues.vcf
  tabix -p vcf -S 1 -f $folder/FuLisF_$data.whole_genome.pvalues.vcf.gz
done

## TajimasD
for data in CEU CHB YRI; do
  gzip -d $folder/TajimasD_$data.whole_genome.pvalues.gz
  cut -d" " -f2- $folder/TajimasD_$data.whole_genome.pvalues > $folder/TajimasD_$data.whole_genome.pvalues_tmp.vcf
  awk -v OFS="\t" '$1=$1' $folder/TajimasD_$data.whole_genome.pvalues_tmp.vcf > $folder/TajimasD_$data.whole_genome.pvalues.vcf
  bgzip $folder/TajimasD_$data.whole_genome.pvalues.vcf
  tabix -p vcf -S 1 -f $folder/TajimasD_$data.whole_genome.pvalues.vcf.gz
done

rm *pvalues
rm *tmp.vcf


### CDTS, bed, 0-based
gzip -d $folder/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed.gz

## recompressing as bgz
bgzip $folder/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed

## indexing
tabix -p bed $folder/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed.gz

