###############################################################################################################################
###############################################################################################################################

### Project: Paralogs 
### Script: Get_paralogs.R 
### Purpose: Get paralogs from Ensembl95 (biomaRt)
### Notes: 1) Only protein coding genes are considered; 2) mapping paralog subtype to  time of duplication event (1:oldest);
### 3) Bidirectional % aa similarity
### Date: 05/03/2019 ##########################################################################################################

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
  left_join(hgnc %>% dplyr::select(hgnc_id,ensembl_gene_id),by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  rename(hgnc_gene_id = hgnc_id) %>%
  left_join(hgnc %>% dplyr::select(hgnc_id,ensembl_gene_id),by = c("hsapiens_paralog_ensembl_gene" = "ensembl_gene_id")) %>%
  rename(hsapiens_paralog_hgnc_gene_id = hgnc_id) %>%
  dplyr::select(ensembl_gene_id,hgnc_gene_id,external_gene_name,hsapiens_paralog_ensembl_gene,
         hsapiens_paralog_hgnc_gene_id,hsapiens_paralog_associated_gene_name,
         hsapiens_paralog_subtype,hsapiens_paralog_orthology_type,hsapiens_paralog_perc_id,
         hsapiens_paralog_perc_id_r1,hsapiens_paralog_paralogy_confidence)



## keep only those paralogues with hgnc id according to hgnc file
## (restrict the associated paralogues to protein coding genes)
## add ordinal phylogeny


subtype <- Cs(Opisthokonta,Bilateria,Chordata,Vertebrata,
              Euteleostomi,Sarcopterygii,Tetrapoda,Amniota,Mammalia,
              Theria,Eutheria,Boreoeutheria,Euarchontoglires,Primates,
              Haplorrhini,Simiiformes,Catarrhini,Hominoidea,Hominidae,
              Homininae,Homo.sapiens)

ensembl95.paralogs.hgnc.proteincoding <- ensembl95.paralogs.hgnc %>%
  mutate(hsapiens_paralog_subtype = ifelse(hsapiens_paralog_subtype =="Homo sapiens","Homo.sapiens",hsapiens_paralog_subtype)) %>%
  mutate(hsapiens_paralog_subtype_ordinal = plyr::mapvalues(hsapiens_paralog_subtype,subtype,c(1:21))) %>%
  dplyr::select(hgnc_gene_id,hsapiens_paralog_hgnc_gene_id,
         hsapiens_paralog_subtype,hsapiens_paralog_subtype_ordinal,
         hsapiens_paralog_orthology_type:hsapiens_paralog_paralogy_confidence) %>%
    dplyr::filter(!is.na(hsapiens_paralog_hgnc_gene_id))


## compute number of paralogues for different % of sequence identity and different
## time of duplication event

paralog.identity.count.unidir <-  function(paralogues,identity) {
  
  count <- paralogues %>% 
    filter(hsapiens_paralog_perc_id >= identity) %>%
    group_by(hgnc_gene_id,hsapiens_paralog_subtype_ordinal) %>% 
    tally() %>% rename(!!paste0("unidir.%identical.aa>",identity) := n)
  return(count)
  
}


paralog.identity.count.bidir <-  function(paralogues,identity) {
  
  count <- paralogues %>% 
    filter(hsapiens_paralog_perc_id >= identity & hsapiens_paralog_perc_id_r1  >= identity) %>%
    group_by(hgnc_gene_id,hsapiens_paralog_subtype_ordinal) %>% 
    tally() %>% rename(!!paste0("bidir.%identical.aa>",identity) := n)
  return(count)
  
}


ensembl95.paralogs.hgnc.proteincoding.count.unidir <- hgnc %>% 
  left_join(paralog.identity.count.unidir(ensembl95.paralogs.hgnc.proteincoding ,0) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)),
    by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count.unidir(ensembl95.paralogs.hgnc.proteincoding ,10) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.unidir(ensembl95.paralogs.hgnc.proteincoding ,20) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.unidir(ensembl95.paralogs.hgnc.proteincoding ,30) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.unidir(ensembl95.paralogs.hgnc.proteincoding ,40) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.unidir(ensembl95.paralogs.hgnc.proteincoding ,50) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::select(-gene_paralog_subtype) %>%
  replace(is.na(.),0)

  
  
ensembl95.paralogs.hgnc.proteincoding.count.bidir <- hgnc %>% 
  left_join(paralog.identity.count.bidir(ensembl95.paralogs.hgnc.proteincoding ,0) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)),
            by=c("hgnc_id"="hgnc_gene_id")) %>%
  left_join(paralog.identity.count.bidir(ensembl95.paralogs.hgnc.proteincoding ,10) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.bidir(ensembl95.paralogs.hgnc.proteincoding ,20) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.bidir(ensembl95.paralogs.hgnc.proteincoding ,30) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.bidir(ensembl95.paralogs.hgnc.proteincoding ,40) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  left_join(paralog.identity.count.bidir(ensembl95.paralogs.hgnc.proteincoding ,50) %>%
              mutate(gene_paralog_subtype = paste0(hgnc_gene_id,"-",hsapiens_paralog_subtype_ordinal)) %>%
              ungroup() %>%
              dplyr::select(4,3),
            by=c("gene_paralog_subtype"="gene_paralog_subtype")) %>%
  dplyr::select(-gene_paralog_subtype) %>%
  replace(is.na(.),0)

################################################################################################################################
################################################################################################################################

## export files

write.table(ensembl95.paralogs.hgnc.proteincoding,
            gzfile("./results/paralogs_ensembl95_proteincoding.txt.gz"),
            quote = F, sep = "\t", row.names = F)

write.table(ensembl95.paralogs.hgnc.proteincoding.count.unidir,
            gzfile("./results/paralogs_count_ensembl95_proteincoding.unidirectional.txt.gz"),
            quote = F, sep = "\t", row.names = F)

write.table(ensembl95.paralogs.hgnc.proteincoding.count.bidir,
            gzfile("./results/paralogs_count_ensembl95_proteincoding.bidirectional.txt.gz"),
            quote = F, sep = "\t", row.names = F)


################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################

