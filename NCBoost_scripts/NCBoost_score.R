# coding: utf-8
# Barthelemy Caron, Clinical BioInformatics Lab, IMAGINE 
# This script will score variants associated to protein coding genes, using the associated gene
# name to choose the corresponding model.

args = commandArgs(trailingOnly=TRUE)
inF <- args[1]
outF <- args[2]
n_rows <- as.integer(args[3])
total_columns <- as.integer(args[4])

options(scipen=999)


library(xgboost)
suppressMessages(library(plyr))
suppressMessages(library(dplyr))

variable_types = list("chr" = "",
                      "pos" = 0, 
                      "ref" = "", 
                      "alt" = "", 
                      "annovar_annotation" = "", 
                      "closest_gene_name"  = "", 
                      "gene_type" = "", 
                      "pLI" = 0, 
                      "familyMemberCount" = 0, 
                      "ncRVIS" = 0, 
                      "ncGERP" = 0, 
                      "RVIS_percentile" = 0, 
                      "slr_dnds" = 0, 
                      "GDI" = 0, 
                      "gene_age" = 0,
                      "UTR3" = 0, 
                      "UTR5" = 0, 
                      "downstream" = 0, 
                      "intergenic" = 0, 
                      "intronic" = 0, 
                      "upstream" = 0, 
                      "GC" = 0, 
                      "CpG"= 0, 
                      "priPhCons" = 0, 
                      "mamPhCons" = 0, 
                      "verPhCons" = 0, 
                      "priPhyloP" = 0, 
                      "mamPhyloP" = 0, 
                      "verPhyloP" = 0, 
                      "GerpN" = 0, 
                      "GerpS" = 0, 
                      "bStatistic" = 0, 
                      "TajimasD_YRI_pvalue" = 0, 
                      "TajimasD_CEU_pvalue" = 0, 
                      "TajimasD_CHB_pvalue" = 0, 
                      "FuLisD_YRI_pvalue" = 0, 
                      "FuLisD_CEU_pvalue" = 0, 
                      "FuLisD_CHB_pvalue" = 0, 
                      "FuLisF_YRI_pvalue"= 0, 
                      "FuLisF_CEU_pvalue" = 0, 
                      "FuLisF_CHB_pvalue" = 0, 
                      "meanDaf1000G" = 0, 
                      "meanHet1000G" = 0, 
                      "meanMAF1000G" = 0, 
                      "meanMAFGnomAD" = 0, 
                      "meanMAF_AFRGnomAD" = 0, 
                      "meanMAF_AMRGnomAD" = 0, 
                      "meanMAF_ASJGnomAD" = 0, 
                      "meanMAF_EASGnomAD" = 0, 
                      "meanMAF_FINGnomAD" = 0, 
                      "meanMAF_NFEGnomAD" = 0, 
                      "meanMAF_OTHGnomAD" = 0, 
                      "CDTS" = 0)

#Loading features
features_selection <- read.table("NCBoost_data/features.tsv", sep="\t", stringsAsFactors = F, header=F, quote="")$V1

#Loading gene2model map
sampled_genes_2_model <- read.table(sprintf("NCBoost_data/gene_to_sample.tsv"), header=T, stringsAsFactors = F, quote="", sep="\t")

#Loading 10-models
for (i in 1:10) {
  e = new.env()
  load(sprintf("NCBoost_data/NCBoost_model_%s.Rdata", i), env=e)
  assign(e$name, e$xgb_model)
  
}

model_names <- c("xgb_model_1","xgb_model_2","xgb_model_3","xgb_model_4","xgb_model_5","xgb_model_6","xgb_model_7","xgb_model_8","xgb_model_9","xgb_model_10")
list_models   <- list(xgb_model_1,xgb_model_2,xgb_model_3,xgb_model_4,xgb_model_5,xgb_model_6,xgb_model_7,xgb_model_8,xgb_model_9,xgb_model_10)

number_of_columns_input <- length(variable_types)

conIN <- file(description=inF, open="r")

cat('\t NCBoost progressing\n\n')
pb <- txtProgressBar(min=1, max=n_rows, style=3)

#scoring & saving
for(i in 1:n_rows) {
  setTxtProgressBar(pb, i)
  if (i == 1) {
    if (total_columns == 53) {
      file_header <- c(scan(file=conIN, nlines=1, quiet=TRUE, sep="\t", quote="", what='char', nmax=number_of_columns_input),
                       'NCBoost')
      standard_names <- file_header[1:number_of_columns_input]
      extra_names <- c()
      input_ncol <- length(file_header)
      list_additional_item_as_char <- vector("list", (total_columns-53)) 

    } else {
      file_header <- c(scan(file=conIN, nlines=1, quiet=TRUE, sep="\t", quote="", what='char', nmax=number_of_columns_input),
                       'NCBoost',
                       scan(file=conIN, nlines=1, quiet=TRUE, sep="\t", quote="", what='char'))
      standard_names <- file_header[1:number_of_columns_input]
      extra_names <- file_header[(number_of_columns_input+3):length(file_header)]
      input_ncol <- length(file_header)
      list_additional_item_as_char <- vector("list", (total_columns-53)) 
      names(list_additional_item_as_char) <- extra_names
      for (item in 1:length(list_additional_item_as_char)) {list_additional_item_as_char[[item]] <- ""}
      
    }
    
    write(matrix(file_header), file=outF, sep="\t", ncolumns=input_ncol)
  } else {
    variant_to_score <- data.frame(scan(file=conIN, nlines=1, quiet=TRUE, sep="\t", quote="", what=append(variable_types, list_additional_item_as_char)), stringsAsFactors = F)
    current_gene <- variant_to_score$closest_gene_name
    associated_sample <- sampled_genes_2_model[sampled_genes_2_model$gene_name == current_gene, "sample_index"]
    
    if (length(associated_sample) != 0) {
      current_features <- variant_to_score[,features_selection]
      variant_to_score[, "NCBoost"] <- predict(list_models[[associated_sample]], as.matrix(current_features), missing=NaN)
      variant_to_score <- variant_to_score[,c(standard_names, "NCBoost", extra_names)]
      write(as.matrix(variant_to_score), file=outF, sep="\t", append=T, ncolumns=input_ncol)
    } else {
      variant_to_score[, "NCBoost"] <- NA
      variant_to_score <- variant_to_score[,c(standard_names, "NCBoost", extra_names)]
      write(as.matrix(variant_to_score), file=outF, sep="\t", append=T, ncolumns=input_ncol)
    }
  }
}
close(pb)
cat('\n \t NCBoost is done.\n')

