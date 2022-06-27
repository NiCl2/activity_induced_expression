# compare Sanchez-Priego activity-induced genes to other gene sets

library(dplyr)
library(biomaRt)

setwd("~/Documents/Projects/activity_induced_expression/")

glu_45m <- read.csv("~/Documents/genesets/ActivityInduced/Sanchez_KCl_Glu-45min.csv") 
glu_45m$GeneSetID = ifelse(glu_45m$log2FoldChange > 0, "glu_45m_up", "glu_45m_down")

glu_4h <- read.csv("~/Documents/genesets/ActivityInduced/Sanchez_KCl_Glu-4h.csv")
glu_4h$GeneSetID = ifelse(glu_4h$log2FoldChange > 0, "glu_4h_up", "glu_4h_down")

gaba_45m <- read.csv("~/Documents/genesets/ActivityInduced/Sanchez_KCl_GABA-45min.csv")
gaba_45m$GeneSetID = ifelse(gaba_45m$log2FoldChange > 0, "gaba_45m_up", "gaba_45m_down")

gaba_4h <- read.csv("~/Documents/genesets/ActivityInduced/Sanchez_KCl_GABA-4h.csv")
gaba_4h$GeneSetID = ifelse(gaba_4h$log2FoldChange > 0, "gaba_4h_up", "gaba_4h_down")

# get EnsemblIDs

ensembl <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
glu_45m_conv <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
              filters='hgnc_symbol', values = glu_45m$Gene, mart = mart)
glu_45m <- left_join(glu_45m, glu_45m_conv, by = c("Gene" = "hgnc_symbol")) %>% 
  filter(!is.na(ensembl_gene_id)) %>% group_by(Gene) %>% filter(n() == 1) %>% ungroup %>% dplyr::select(GeneSetID, ensembl_gene_id)

glu_4h_conv <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                      filters='hgnc_symbol', values = glu_4h$Gene, mart = mart)
glu_4h <- left_join(glu_4h, glu_4h_conv, by = c("Gene" = "hgnc_symbol")) %>% 
  filter(!is.na(ensembl_gene_id)) %>% group_by(Gene) %>% filter(n() == 1) %>% ungroup %>% dplyr::select(GeneSetID, ensembl_gene_id)

gaba_45m_conv <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                      filters='hgnc_symbol', values = gaba_45m$Gene, mart = mart)
gaba_45m <- left_join(gaba_45m, gaba_45m_conv, by = c("Gene" = "hgnc_symbol")) %>% 
  filter(!is.na(ensembl_gene_id)) %>% group_by(Gene) %>% filter(n() == 1) %>% ungroup %>% dplyr::select(GeneSetID, ensembl_gene_id)

gaba_4h_conv <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                      filters='hgnc_symbol', values = gaba_4h$Gene, mart = mart)
gaba_4h <- left_join(gaba_4h, gaba_4h_conv, by = c("Gene" = "hgnc_symbol")) %>% 
  filter(!is.na(ensembl_gene_id)) %>% group_by(Gene) %>% filter(n() == 1) %>% ungroup %>% dplyr::select(GeneSetID, ensembl_gene_id) 

# Roussos KCl-induced gene expression (old and newer data)
roussos_kcl1 <- read.table("~/Documents/genesets/ActivityInduced/Roussos_KCl_Activity-induced.txt", col.names = c("GeneSetID", "ensembl_gene_id"))
roussos_kcl2 <- read.table("~/Documents/genesets/ActivityInduced/KCl_genesets.txt", col.names = c("GeneSetID", "ensembl_gene_id")) %>%
  filter(GeneSetID == "up1" | GeneSetID == "up6" | GeneSetID == "down1" | GeneSetID == "down6")

# LTP-induced expression and TRAP
ltp <- read.table("~/Documents/RNAseq/LTP/LTP_ensembl.txt", col.names = c("GeneSetID", "ensembl_gene_id"))

# FMRP targets
fmrp <- read.table("~/Documents/genesets/FMRP/FMRPtargetsDarnellFDR001_ensembl.txt", col.names = c("GeneSetID", "ensembl_gene_id"))



genesets <- rbind(glu_45m, glu_4h, gaba_45m, gaba_4h, roussos_kcl1, roussos_kcl2, ltp, fmrp)

# intersect all
intersect_all <- function(df) {
  gs_matrix = data.frame(row.names = unique(df$GeneSetID))
  for (gs1 in unique(df$GeneSetID)) {
    for (gs2 in unique(df$GeneSetID)) {
      gs_matrix[gs2, gs1] <- length(intersect(df$ensembl_gene_id[df$GeneSetID == gs1], df$ensembl_gene_id[df$GeneSetID == gs2]))
    }
  }
  return(gs_matrix)
}



