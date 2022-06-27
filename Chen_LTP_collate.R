# collate LTP data from Chen et al 2017

setwd("~/Documents/RNAseq/LTP")

library(dplyr)

# 30 minutes
ltp30 <- read.csv("GSE79790_30_min_data.csv") %>% 
  dplyr::select(X, entrez, logFC_TRAP30LTP_vs_TRAP30basal, FDR_TRAP30LTP_vs_TRAP30basal, logFC_Total30LTP_vs_Total30basal, FDR_Total30LTP_vs_Total30basal) %>%
  mutate(total_up = FDR_Total30LTP_vs_Total30basal < 0.05 & logFC_Total30LTP_vs_Total30basal > 0) %>% 
  mutate(total_down = FDR_Total30LTP_vs_Total30basal < 0.05 & logFC_Total30LTP_vs_Total30basal < 0) %>%
  mutate(trap_up = FDR_TRAP30LTP_vs_TRAP30basal < 0.05 & logFC_TRAP30LTP_vs_TRAP30basal > 0) %>% 
  mutate(trap_down = FDR_TRAP30LTP_vs_TRAP30basal < 0.05 & logFC_TRAP30LTP_vs_TRAP30basal < 0)
ltp30_gs <- data.frame(GeneSetID = "ltp30_total_up", EntrezID = ltp30$entrez[ltp30$total_up %in% TRUE])
ltp30_gs <- rbind(ltp30_gs, 
                  data.frame(GeneSetID = "ltp30_total_down", EntrezID = ltp30$entrez[ltp30$total_down %in% TRUE]))
ltp30_gs <- rbind(ltp30_gs, 
                  data.frame(GeneSetID = "ltp30_trap_up", EntrezID = ltp30$entrez[ltp30$trap_up %in% TRUE]))
ltp30_gs <- rbind(ltp30_gs, 
                  data.frame(GeneSetID = "ltp30_trap_down", EntrezID = ltp30$entrez[ltp30$trap_down %in% TRUE]))

# 60 minutes
ltp60 <- read.csv("GSE79790_60_min_data.csv") %>%
  dplyr::select(X, entrez, logFC_TRAP60LTP_vs_TRAP60basal, FDR_TRAP60LTP_vs_TRAP60basal, logFC_Total60LTP_vs_Total60basal, FDR_Total60LTP_vs_Total60basal) %>%
  mutate(total_up = FDR_Total60LTP_vs_Total60basal < 0.05 & logFC_Total60LTP_vs_Total60basal > 0) %>% 
  mutate(total_down = FDR_Total60LTP_vs_Total60basal < 0.05 & logFC_Total60LTP_vs_Total60basal < 0) %>%
  mutate(trap_up = FDR_TRAP60LTP_vs_TRAP60basal < 0.05 & logFC_TRAP60LTP_vs_TRAP60basal > 0) %>% 
  mutate(trap_down = FDR_TRAP60LTP_vs_TRAP60basal < 0.05 & logFC_TRAP60LTP_vs_TRAP60basal < 0)
ltp60_gs <- data.frame(GeneSetID = "ltp60_total_up", EntrezID = ltp60$entrez[ltp60$total_up %in% TRUE])
ltp60_gs <- rbind(ltp60_gs, 
                  data.frame(GeneSetID = "ltp60_total_down", EntrezID = ltp60$entrez[ltp60$total_down %in% TRUE]))
ltp60_gs <- rbind(ltp60_gs, 
                  data.frame(GeneSetID = "ltp60_trap_up", EntrezID = ltp60$entrez[ltp60$trap_up %in% TRUE]))
ltp60_gs <- rbind(ltp60_gs, 
                  data.frame(GeneSetID = "ltp60_trap_down", EntrezID = ltp60$entrez[ltp60$trap_down %in% TRUE]))

# 120 minutes
ltp120 <- read.csv("GSE79790_120_min_data.csv") %>%
  dplyr::select(X, entrez, logFC_TRAP120LTP_vs_TRAP120basal, FDR_TRAP120LTP_vs_TRAP120basal, logFC_Total120LTP_vs_Total120basal, FDR_Total120LTP_vs_Total120basal) %>%
  mutate(total_up = FDR_Total120LTP_vs_Total120basal < 0.05 & logFC_Total120LTP_vs_Total120basal > 0) %>% 
  mutate(total_down = FDR_Total120LTP_vs_Total120basal < 0.05 & logFC_Total120LTP_vs_Total120basal < 0) %>%
  mutate(trap_up = FDR_TRAP120LTP_vs_TRAP120basal < 0.05 & logFC_TRAP120LTP_vs_TRAP120basal > 0) %>% 
  mutate(trap_down = FDR_TRAP120LTP_vs_TRAP120basal < 0.05 & logFC_TRAP120LTP_vs_TRAP120basal < 0)
ltp120_gs <- data.frame(GeneSetID = "ltp120_total_up", EntrezID = ltp120$entrez[ltp120$total_up %in% TRUE])
ltp120_gs <- rbind(ltp120_gs, 
                  data.frame(GeneSetID = "ltp120_total_down", EntrezID = ltp120$entrez[ltp120$total_down %in% TRUE]))
ltp120_gs <- rbind(ltp120_gs, 
                  data.frame(GeneSetID = "ltp120_trap_up", EntrezID = ltp120$entrez[ltp120$trap_up %in% TRUE]))
ltp120_gs <- rbind(ltp120_gs, 
                  data.frame(GeneSetID = "ltp120_trap_down", EntrezID = ltp120$entrez[ltp120$trap_down %in% TRUE]))

ltp_all_gs <- rbind(ltp30_gs, ltp60_gs, ltp120_gs) %>% filter(!is.na(EntrezID))

library(biomaRt)

ensembl <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
glu_45m_conv <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                      filters='hgnc_symbol', values = glu_45m$Gene, mart = mart)