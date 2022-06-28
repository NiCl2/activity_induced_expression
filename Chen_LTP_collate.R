# collate LTP data from Chen et al 2017

setwd("~/Documents/RNAseq/LTP")

library(dplyr)

# 30 minutes
# Differential expression cut off defined in paper: P<0.1 & log2fc > 0.4
# Differential expression cut off defined in Meeson thesis: P<0.01
ltp30 <- read.csv("GSE79790_30_min_data.csv") %>% 
  dplyr::select(X, entrez, logFC_TRAP30LTP_vs_TRAP30basal, FDR_TRAP30LTP_vs_TRAP30basal, logFC_Total30LTP_vs_Total30basal, FDR_Total30LTP_vs_Total30basal) %>%
  mutate(total_up = FDR_Total30LTP_vs_Total30basal < 0.01 & logFC_Total30LTP_vs_Total30basal > 0) %>% 
  mutate(total_down = FDR_Total30LTP_vs_Total30basal < 0.01 & logFC_Total30LTP_vs_Total30basal < 0) %>%
  mutate(trap_up = FDR_TRAP30LTP_vs_TRAP30basal < 0.01 & logFC_TRAP30LTP_vs_TRAP30basal > 0) %>% 
  mutate(trap_down = FDR_TRAP30LTP_vs_TRAP30basal < 0.01 & logFC_TRAP30LTP_vs_TRAP30basal < 0)
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
  mutate(total_up = FDR_Total60LTP_vs_Total60basal < 0.01 & logFC_Total60LTP_vs_Total60basal > 0) %>% 
  mutate(total_down = FDR_Total60LTP_vs_Total60basal < 0.01 & logFC_Total60LTP_vs_Total60basal < 0) %>%
  mutate(trap_up = FDR_TRAP60LTP_vs_TRAP60basal < 0.01 & logFC_TRAP60LTP_vs_TRAP60basal > 0) %>% 
  mutate(trap_down = FDR_TRAP60LTP_vs_TRAP60basal < 0.01 & logFC_TRAP60LTP_vs_TRAP60basal < 0)
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
  mutate(total_up = FDR_Total120LTP_vs_Total120basal < 0.01 & logFC_Total120LTP_vs_Total120basal > 0) %>% 
  mutate(total_down = FDR_Total120LTP_vs_Total120basal < 0.01 & logFC_Total120LTP_vs_Total120basal < 0) %>%
  mutate(trap_up = FDR_TRAP120LTP_vs_TRAP120basal < 0.01 & logFC_TRAP120LTP_vs_TRAP120basal > 0) %>% 
  mutate(trap_down = FDR_TRAP120LTP_vs_TRAP120basal < 0.01 & logFC_TRAP120LTP_vs_TRAP120basal < 0)
ltp120_gs <- data.frame(GeneSetID = "ltp120_total_up", EntrezID = ltp120$entrez[ltp120$total_up %in% TRUE])
ltp120_gs <- rbind(ltp120_gs, 
                  data.frame(GeneSetID = "ltp120_total_down", EntrezID = ltp120$entrez[ltp120$total_down %in% TRUE]))
ltp120_gs <- rbind(ltp120_gs, 
                  data.frame(GeneSetID = "ltp120_trap_up", EntrezID = ltp120$entrez[ltp120$trap_up %in% TRUE]))
ltp120_gs <- rbind(ltp120_gs, 
                  data.frame(GeneSetID = "ltp120_trap_down", EntrezID = ltp120$entrez[ltp120$trap_down %in% TRUE]))

# background
bkg.30 <- read.csv("GSE79790_30_min_data.csv") %>% dplyr::select(entrez, Total30basal1_normalizedCount, Total30basal2_normalizedCount, Total30basal3_normalizedCount) %>% filter(!is.na(entrez))
bkg.60 <- read.csv("GSE79790_60_min_data.csv") %>% dplyr::select(entrez, Total60basal1_normalizedCount, Total60basal2_normalizedCount, Total60basal3_normalizedCount) %>% filter(!is.na(entrez))
bkg.120 <- read.csv("GSE79790_120_min_data.csv") %>% dplyr::select(entrez, Total120basal1_normalizedCount, Total120basal2_normalizedCount, Total120basal3_normalizedCount) %>% filter(!is.na(entrez))
bkg <- left_join(bkg.30, bkg.60, by = "entrez") %>% left_join(bkg.120, by = "entrez") 
bkg_filter_index <- rowSums(bkg[, 2:10] > 10) >= 3 # genes where the normalised count is greater than 10 in 3 or more basal samples
bkg <- filter(bkg, bkg_filter_index)
bkg <- data.frame(GeneSetID = "ltp_background", EntrezID = bkg$entrez)

# merge
ltp_all_gs <- rbind(ltp30_gs, ltp60_gs, ltp120_gs, bkg) %>% filter(!is.na(EntrezID))

# convert to human ensembl IDs
library(biomaRt)
mm <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
hu <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
gs_conv_hu <- getLDS(attributes = c("entrezgene_id"), 
                     filters="entrezgene_id", values = ltp_all_gs$EntrezID, mart = mm, 
                     attributesL = c("ensembl_gene_id"), 
                     martL = hu)

ltp_all_gs <- left_join(ltp_all_gs, gs_conv_hu, by = c("EntrezID" = "NCBI.gene..formerly.Entrezgene..ID"))
ltp_all_gs <- filter(ltp_all_gs, !is.na(Gene.stable.ID)) %>% 
  group_by(GeneSetID, EntrezID) %>% filter(n() == 1) %>% ungroup %>%
  dplyr::select(GeneSetID, Gene.stable.ID)

write.table(ltp_all_gs, "Chen_LTP_hu_ensembl.txt", row.names = F, col.names = F, quote = F, sep = "\t")

