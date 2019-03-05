###############################################################################################################################
###############################################################################################################################

### Project: Paralogs #########################################################################################################
### Script: Get_paralogs.R ####################################################################################################
### Purpose: Get paralogs from Ensembl95 (biomaRt) ############################################################################
### Author: Pilar Cacheiro ####################################################################################################
### Date: 15/03/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr)

################################################################################################################################
################################################################################################################################

## import files ################################################################################################################

## hgnc protein coding genes

hgnc <-  read_delim("./data/gene_with_protein_product.txt",delim="\t") %>%
  select(hgnc_id,symbol,ensembl_gene_id)


################################################################################################################################
################################################################################################################################

## retrieve ensembl paralogs - biomart ensembl genes 95 -

library(biomaRt)

mart.human.ensembl95 <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")



ensembl95.paralogs <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                           "hsapiens_paralog_ensembl_gene",
                                           "hsapiens_paralog_associated_gene_name",
                                           "hsapiens_paralog_chromosome",
                                           "hsapiens_paralog_subtype",
                                           "hsapiens_paralog_orthology_type",
                                           "hsapiens_paralog_perc_id",
                                           "hsapiens_paralog_perc_id_r1",
                                           "hsapiens_paralog_paralogy_confidence"),
                            filters = "ensembl_gene_id",
                            values = hgnc$ensembl_gene_id,
                            mart = mart.human.ensembl95)


## merge with hgnc file to get hgnc id 

ensembl95.paralogs.hgnc <-  ensembl95.paralogs %>%
  left_join(hgnc %>% select(hgnc_id,ensembl_gene_id),by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  rename(hgnc_gene_id = hgnc_id) %>%
  left_join(hgnc %>% select(hgnc_id,ensembl_gene_id),by = c("hsapiens_paralog_ensembl_gene" = "ensembl_gene_id")) %>%
  rename(hsapiens_paralog_hgnc_gene_id = hgnc_id) %>%
  select(ensembl_gene_id,hgnc_gene_id,external_gene_name,hsapiens_paralog_ensembl_gene,
         hsapiens_paralog_hgnc_gene_id,hsapiens_paralog_associated_gene_name,
         hsapiens_paralog_subtype,hsapiens_paralog_orthology_type,hsapiens_paralog_perc_id,
         hsapiens_paralog_perc_id_r1,hsapiens_paralog_paralogy_confidence)



## keep only those paralogues with hgnc id according to hgnc file
## (restrict the associated paralogues to protein coding genes)


ensembl95.paralogs.hgnc.proteincoding <- ensembl95.paralogs.hgnc %>%
  select(hgnc_gene_id,hsapiens_paralog_hgnc_gene_id,
         hsapiens_paralog_subtype:hsapiens_paralog_paralogy_confidence) %>%
  filter(!is.na(hsapiens_paralog_hgnc_gene_id))


## compute number of paralogues for different % of sequence identity

paralog.identity.count <-  function(paralogues,identity) {
  
  count <- paralogues %>% 
    filter(hsapiens_paralog_perc_id >= identity) %>%
    group_by(hgnc_gene_id) %>% 
    tally() %>% rename(!!paste0("% identical aa >",identity) := n)
  return(count)
  
}

ensembl95.paralogs.hgnc.proteincoding.count <- hgnc %>% 
  left_join(paralog.identity.count(ensembl95.paralogs.hgnc.proteincoding ,0),by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count(ensembl95.paralogs.hgnc.proteincoding ,10),by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count(ensembl95.paralogs.hgnc.proteincoding ,20),by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count(ensembl95.paralogs.hgnc.proteincoding ,30),by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count(ensembl95.paralogs.hgnc.proteincoding ,40),by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count(ensembl95.paralogs.hgnc.proteincoding ,50),by=c("hgnc_id"="hgnc_gene_id")) %>%
  replace(is.na(.),0)

################################################################################################################################
################################################################################################################################

## export files

write.table(ensembl95.paralogs.hgnc.proteincoding,"./results/paralogs_ensembl95_proteincoding.txt",
            quote = F, sep = "\t", row.names = F)

write.table(ensembl95.paralogs.hgnc.proteincoding.count,"./results/paralogs_count_ensembl95_proteincoding.txt",
            quote = F, sep = "\t", row.names = F)

################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

