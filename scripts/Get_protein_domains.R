###############################################################################################################################
###############################################################################################################################

### Project: Paralogs 
### Script: Get_protein_domains.R 
### Purpose: Get interpro and pfam annotations (functional domains)
### Notes: Retrieve information on functional domains (first: interpro and pfam) for all protein coding genes
### Date: 07/03/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr);library(Hmisc)

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################

## hgnc protein coding genes

hgnc <-  read_delim("./data/gene_with_protein_product.txt",delim="\t") %>%
  select(hgnc_id,symbol,ensembl_gene_id)


################################################################################################################################
################################################################################################################################

## retrieve functional domains - biomart ensembl genes 95 -

library(biomaRt)

mart.human.ensembl95 <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")


ensembl95.interpro <- getBM(attributes = c("hgnc_id","interpro"),
                            filters = "hgnc_id",
                            values = hgnc$hgnc_id,
                            mart = mart.human.ensembl95) %>%
  filter(interpro!="") %>%
  distinct()


ensembl95.pfam <- getBM(attributes = c("hgnc_id","pfam"),
                            filters = "hgnc_id",
                            values = hgnc$hgnc_id,
                            mart = mart.human.ensembl95) %>%
  filter(pfam!="") %>%
  distinct()

## many other domain databases : pirsf, prints, scanprosite, pfscan, smart, superfamily, tigrfam

################################################################################################################################
################################################################################################################################

## export files


write.table(ensembl95.interpro ,
            gzfile("./results/ensembl95.interpro.txt.gz"),
            quote = F, sep = "\t", row.names = F)

write.table(ensembl95.pfam,
            gzfile("./results/ensembl95.pfam.txt.gz"),
            quote = F, sep = "\t", row.names = F)
                                          