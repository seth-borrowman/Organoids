library(tximeta)
library(DESeq2)
library(tidyverse)
library(ashr)
library(biomaRt)
library(ggrepel)
library(tximport)

# Set wd
setwd("Z:/Organoids/RNAseq/quants")
#ifelse(!dir.exists('DE'), dir.create('DE'), FALSE)

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

files_needed <- c("Alz_iPSC_quant", "Alz_Org_quant",
                  "HFF_iPSC_quant", "HFF_Org_quant")
condition <- c("Case", "Case", "Control", "Control")

counts <- getBM(attributes = "hgnc_symbol", mart = mart)

for (i in files_needed) {
    # Load case data
    case <- tximport(paste(i, "/quant.sf", sep = ""),
                     type = "salmon",
                     txOut = T) %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        dplyr::rename(ensembl_transcript_id_version = rowname,
                      NumReads = counts)
    
    geneNames <- getBM(filters = "ensembl_transcript_id_version",
                       attributes = c("ensembl_transcript_id_version",
                                      "ensembl_gene_id_version",
                                      "hgnc_symbol"),
                       values = case$ensembl_transcript_id_version,
                       mart = mart)
    case <- merge(case, geneNames)
    
    # IIHG gives results by transcript ID, want by gene - convert to gene
    case_sum <- case %>%
        group_by(hgnc_symbol) %>%
        summarise(Reads = sum(NumReads))
    colnames(case_sum)[2] <- i
    
    counts <- merge(counts, case_sum, by = "hgnc_symbol")
}

row.names(counts) <- counts[,1]
counts <- counts[,-1]

ddsTxi <- DESeqDataSetFromMatrix(round(counts), DataFrame(condition),
                                 ~ condition, tidy = F)

# DESeq analysis with ASHR shrinkage
dds <- DESeq(ddsTxi)
resLFC <- lfcShrink(dds, contrast = c("condition", "Case", "Control"), type = "ashr")
resultsDF <- as.data.frame(resLFC)
resultsDF$GeneID <- resLFC@rownames

resultsDF <- resultsDF %>%
    mutate(sig = case_when(
        -log(padj) >= -log(0.05/nrow(resultsDF)) &
            abs(log2FoldChange) >= 2 ~ 1,
        GeneID == "HIPK2" ~ 1,
        .default = 0
    ))
resultsDF$sig <- factor(resultsDF$sig, ordered = is.ordered(c(0, 1)))

options(ggrepel.max.overlaps = Inf)
ggplot(data = resultsDF, aes(x = log2FoldChange, y = -log(padj), color = sig)) +
    geom_point(show.legend = F) +
    #geom_vline(xintercept = c(-2, 2), color = "pink") +
    geom_hline(yintercept = -log(0.05/nrow(resultsDF)), color = "lightblue") +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    geom_text_repel(data = resultsDF %>%
                        filter(sig == 1),
                    aes(label = GeneID)) +
    theme(legend.position = "none") +
    ggtitle("Alz vs HFF")

### ----
condition <- c("Control", "Case", "Control", "Case")

ddsTxi <- DESeqDataSetFromMatrix(round(counts), DataFrame(condition),
                                 ~ condition, tidy = F)

# DESeq analysis with ASHR shrinkage
dds <- DESeq(ddsTxi)
resLFC <- lfcShrink(dds, contrast = c("condition", "Case", "Control"), type = "ashr")
resultsDF <- as.data.frame(resLFC)
resultsDF$GeneID <- resLFC@rownames

resultsDF <- resultsDF %>%
    mutate(sig = case_when(
        -log(padj) >= -log(0.05/nrow(resultsDF)) &
            abs(log2FoldChange) >= 2 ~ 1,
        .default = 0
    ))
resultsDF$sig <- factor(resultsDF$sig, ordered = is.ordered(c(0, 1)))

options(ggrepel.max.overlaps = 10)
ggplot(data = resultsDF, aes(x = log2FoldChange, y = -log(padj), color = sig)) +
    geom_point(show.legend = F) +
    #geom_vline(xintercept = c(-2, 2), color = "pink") +
    geom_hline(yintercept = -log(0.05/nrow(resultsDF)), color = "lightblue") +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    geom_text_repel(data = resultsDF %>%
                        filter(sig == 1),
                    aes(label = GeneID)) +
    theme(legend.position = "none") +
    ggtitle("Organoids vs iPSCs")
