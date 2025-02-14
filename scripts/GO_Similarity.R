##############################################################################################################################

### Project: Paralogs 
### Script: GO_Similiarity.R 

###############################################################################################################################
###############################################################################################################################

## load packages ##############################################################################################################

library(dplyr);library(tidyr);library(stringr);library(readr)
library(GOSemSim)

################################################################################################################################
################################################################################################################################


## import files ################################################################################################################

gene.identifiers <- read_delim("./data/gene_with_protein_product.txt",delim="\t") %>%
  dplyr::select(hgnc_id,entrez_id)

paralogs <- read_delim("./results/paralogs_ensembl95_proteincoding.txt.gz",delim="\t") %>%
  inner_join(gene.identifiers,by=c("hgnc_gene_id" = "hgnc_id")) %>%
  rename(entrez_gene = entrez_id) %>%
  inner_join(gene.identifiers,by=c("hsapiens_paralog_hgnc_gene_id" = "hgnc_id")) %>%
  rename(entrez_paralog = entrez_id)

paralog.pairs <- paralogs %>%
  dplyr::select(entrez_gene,entrez_paralog) %>%
  distinct()


## import GO annotations ####################################################################################################

hsbp <- godata('org.Hs.eg.db', ont="BP")
hsmf <- godata('org.Hs.eg.db', ont="MF")
hscc <- godata('org.Hs.eg.db', ont="CC")
hsbp.go <- hsbp@geneAnno

pair.go.similarity.df <- list()

for(i in 1:dim(paralog.pairs)[1]){
  
  for(i in 1:10){
    
    gene = paralog.pairs[i,1]
    paralog = paralog.pairs[i,2]
    
    sim.bp.reisnik = mgeneSim(genes = c(gene,paralog), semData =hsbp, measure = "Resnik", drop = "IEA", combine = "BMA",verbose = TRUE) 
    sim.mf.reisnik = mgeneSim(genes = c(gene,paralog), semData =hsmf, measure = "Resnik", drop = "IEA", combine = "BMA",verbose = TRUE) 
    sim.cc.reisnik = mgeneSim(genes = c(gene,paralog), semData =hscc, measure = "Resnik", drop = "IEA", combine = "BMA",verbose = TRUE) 
    
    if(length(sim.bp.reisnik)>1) {sim.bp.reisnik.score = sim.bp.reisnik[1,2] 
    } else {
      sim.bp.reisnik.score = NA
    }
    
    if(length(sim.mf.reisnik)>1) {sim.mf.reisnik.score = sim.mf.reisnik[1,2] 
    } else {
      sim.mf.reisnik.score = NA
    }
    
    if(length(sim.cc.reisnik)>1) {sim.cc.reisnik.score = sim.cc.reisnik[1,2] 
    } else {
      sim.cc.reisnik.score = NA
    }
    
    
    pair.go.similarity.df[[i]] = data.frame(gene,paralog,sim.bp.reisnik.score,sim.mf.reisnik.score,sim.cc.reisnik.score)
    
    cat("iteration:", i, " \n") 
    
    flush.console()
    
  }
  
pair.go.similarity.df.all <- do.call(rbind,pair.go.similarity.df)