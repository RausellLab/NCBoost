# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 

# create a folder to gather the downloaded data
storage_folder='NCBoost_features'
mkdir $storage_folder

### Download CADD v1.3 prescored whole genome 
wget --directory-prefix=$storage_folder http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz
wget --directory-prefix=$storage_folder http://krishna.gs.washington.edu/download/CADD/v1.3/whole_genome_SNVs_inclAnno.tsv.gz.tbi

### Download 1000 genome Project Tajima's D, FuLi's F* and D*
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisD_CEU.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisD_CHB.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisD_YRI.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisF_CEU.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisF_CHB.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/FuLisF_YRI.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/TajimasD_CEU.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/TajimasD_CHB.whole_genome.pvalues.gz
wget --directory-prefix=$storage_folder http://hsb.upf.edu/hsb_data/positive_selection_NAR2013/TajimasD_YRI.whole_genome.pvalues.gz

### Download CDTS score
wget --directory-prefix=$storage_folder http://www.hli-opendata.com/noncoding/Pipeline/CDTS_diff_perc_coordsorted_gnomAD_N15496_hg19.bed.gz
 
### Download mean DAF and HET from GWAVA
wget --directory-prefix=$storage_folder ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/source_data/1kg/daf.bed.gz
wget --directory-prefix=$storage_folder ftp://ftp.sanger.ac.uk/pub/resources/software/gwava/v1.0/source_data/1kg/het_rates.bed.gz

### Download GnomAD indexed and corresponding index files
for chr in {1..22} X; do
	wget --directory-prefix=$storage_folder https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr$chr.vcf.bgz
	wget --directory-prefix=$storage_folder https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr$chr.vcf.bgz.tbi
done
