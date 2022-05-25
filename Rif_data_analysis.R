## Analyze MeDIPseq data: ---------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("qsea")
#BiocManager::install("BSgenome")
library("BSgenome")
#available.genomes()
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(qsea)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(tidyverse)
library(annotatr)
source("./data/QSEA_outcome_hepatic_Rif/R_functions.R")
load("./data/QSEA_outcome_hepatic_Rif/ROIs_2_2021Jan19.RData")
load("./data/QSEA_outcome_hepatic_Rif/annotations_2021Sep24.RData")

# Extract and annotate DMRs for DPI samples in each dose (Therapeutic & Toxix)
load("./data/QSEA_outcome_hepatic_Rif/qsea_outcome_Rif_per_dose_20220407.RData")
names(outcome)

MeDIP_Rif_dose <- list()
for (treat in names(outcome)) {
  sig <- isSignificant(outcome[[treat]]$qseaGLM, fdr_th = 0.01)
  QSEA_outcome <- makeTable(outcome[[treat]]$qseaSet, glm = outcome[[treat]]$qseaGLM,
                            groupMeans=getSampleGroups(outcome[[treat]]$qseaSet),
                            keep=sig, annotation = c(ROIs_2), norm_methods = "beta")
  annotated_genes <- get.annotated_genes(QSEA_outcome, annotations)
  
  gene_DMRs <-  annotated_genes[which(annotated_genes$DMR_count>quantile(annotated_genes$DMR_count, .95)),]
  genes_DMRs_promoter <- gene_DMRs[which(gene_DMRs$DMR_in_promoter_count>quantile(gene_DMRs$DMR_in_promoter_count, .95)),] 
  avg_genes_DMRs_protomer <- get.avg_DMRs_genes(genes_DMRs_promoter, QSEA_outcome)
  selected_genes <- avg_genes_DMRs_protomer[which(abs(avg_genes_DMRs_protomer$log2FC_avg)>0.5),]
  write.csv(selected_genes, file = paste0("MeDIP_selectedGenes_", treat, "_20220407.csv"))
  
  avg_genes <- get.avg_DMRs_genes(annotated_genes, QSEA_outcome)
  volcano_gene <- pre.volcano_plot(avg_genes, selected_genes$annot.symbol, Log2FC_cutoff =0.5)
  
  ggplot(annotated_genes, aes(x=DMR_count)) + geom_histogram(binwidth = 1, color="gray") + 
    geom_histogram(data=gene_DMRs, aes(x=DMR_count), binwidth=1, color="transparent", fill="red")
  ggplot(gene_DMRs, aes(x=DMR_in_promoter_count)) + geom_histogram(binwidth = 1) +
    geom_histogram(data=genes_DMRs_promoter, aes(x=DMR_in_promoter_count), binwidth = 1, color="transparent", fill="red")
  
  MeDIP_Rif_dose[[treat]] <- list("sig" = sig, "QSEA_outcome" = QSEA_outcome, "annotated_genes" = annotated_genes,
                            "gene_DMRs" = gene_DMRs, "genes_DMRs_promoter" = genes_DMRs_promoter, "avg_genes_DMRs_protomer" = avg_genes_DMRs_protomer, 
                            "selected_genes" = selected_genes, "avg_genes" = avg_genes, "volcano_gene" = volcano_gene)  
  rm(sig, QSEA_outcome, annotated_genes, gene_DMRs, genes_DMRs_promoter,
     avg_genes_DMRs_protomer, selected_genes, avg_genes, volcano_gene)
}
names(MeDIP_Rif_dose)

# make volcano plots
pdf(paste0("./outcome/MeDIP_Rif_doses_20220408.pdf"), onefile=TRUE)
get.volcano_plot(MeDIP_Rif_dose$Rif_The$volcano_gene, Log2FC_cutoff = 0.5)
get.volcano_plot(MeDIP_Rif_dose$Rif_Tox$volcano_gene, Log2FC_cutoff = 0.5)
dev.off()

## Analyze RNAseq data: ---------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggfortify)

load("./data/RNAseq_data/Rif_RNAseq_DEseq2_allsamples.RData")

#overview of the RNAseq - PCA plot
pdf("./outcome/RNAseq_pca.pdf", height = 4, width =6)
autoplot(DESeq2::plotPCA(norm_data, intgroup = "condition", ntop = 40000)) + theme_bw()
dev.off()

# tidy the RNAseq normalized data:
rna_annot_info <- read.table(file = "./data/RNAseq_data/sample_info.txt",
                             sep = "\t", header = TRUE, colClasses = 'character') %>%
  mutate(sample_name_id = paste0(condition, "_", sampling_time_point, "_", sample_ID))

RNAseq_norm_dat <- assay(norm_data) %>% t() %>% as.data.frame() %>%
  rownames_to_column(var = "sample_name_replicate") %>%
  left_join(rna_annot_info[, c("sample_name_replicate", "sample_name_id")]) %>%
  select(-c("sample_name_replicate")) %>% 
  column_to_rownames("sample_name_id") %>% t() %>% as.data.frame() %>%
  rename_with(~gsub("ConDMSO", "Control", .x))# ajust the sample names


## Analyze protein data: ---------------------------------------------
library(tidyverse)
# load and tidy the data
input_dat <- read.csv("./data/protein_data/Hecatos_Hepatic_Px_RIF_Ther_Tox_log2_norm.txt", 
                      sep = "\t") %>% # get Rifampicin (RIF) normalized data
  full_join(read.csv("./data/protein_data/Hecatos_Hepatic_Px_01DMSO_log2_normalized.txt",
                     sep = "\t")) %>% # get control normalized data
  filter(., !grepl(":" , Row.Names)) %>% # only select protein mapped to unique ID
  column_to_rownames("Row.Names") %>%
  filter(!is.na(rowSums(., na.rm = TRUE))) # remove protein with NA in all samples

annot_info <- read.csv("./data/protein_data/Hecatos_Hepatic_Px_RIF_Ther_Tox_SampleAnnot.txt",
                       sep = "\t") %>% # get sample annotation information
  select("Fraction", "Treatment...Time.Point...ID")
annot_info$Name <- paste0("X", annot_info$Fraction)

protein_norm_dat <- as.data.frame(t(input_dat)) %>% rownames_to_column(var = "Name") %>%
  left_join(., annot_info, by="Name") %>% select(!c("Fraction", "Name")) %>% # convert fraction IDs to sample names
  column_to_rownames( "Treatment...Time.Point...ID") %>% 
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "Annot_info") %>% # extract the Uniprot IDs from protein annotation step
  mutate("UNIPROT" = str_split_fixed(Annot_info, "[|]", 3)[,2]) %>%
  select(!c("Annot_info")) %>% column_to_rownames("UNIPROT") %>%
  rename_with(~gsub("01%DMSO", "Control", .x)) %>% # ajust the sample names - treatment conditions
  rename_with(~gsub("RIF_Ther", "Rif_The", .x)) %>%
  rename_with(~gsub("RIF_Tox", "Rif_Tox", .x)) %>%
  rename_with(~gsub("T336", "336", .x)) %>% # ajust the sample names - time points
  rename_with(~gsub("T240", "240", .x)) %>%
  rename_with(~gsub("T168", "168", .x)) %>%
  rename_with(~gsub("T72", "072", .x)) %>%
  rename_with(~gsub("T24", "024", .x)) %>%
  rename_with(~gsub("T8", "008", .x)) %>%
  rename_with(~gsub("T2", "002", .x))
# overview of the proteomics data - PCA plot
library(factoextra)
protein_pca <- as.data.frame(t(protein_norm_dat))
protein_pca[is.na(protein_pca)] <- 0
res.pca <- prcomp(protein_pca)
protein_pca$condition <- substr(rownames(protein_pca), 1, 7)

pdf("./outcome/protein_pca.pdf", height = 3, width = 6)
autoplot(res.pca, data=protein_pca, 
         colour = "condition", size = 3)+ theme_bw()
dev.off()

## compare clustering profile: RNAseq vs. proteomics data --------------------------------------
library(dendextend)
library("scales")
dend_RNAseq <- t(RNAseq_norm_dat) %>% scale %>% dist %>% 
  hclust("ward.D") %>% as.dendrogram

dend_protein <- t(protein_norm_dat) %>% scale %>% dist %>% 
  hclust("ward.D") %>% as.dendrogram

groupCodes <- substr(colnames(protein_norm_dat), 1, 7)
colorCodes <- c( "#CD9600", "#C77CFF", "#d87876")
names(colorCodes) <- unique(groupCodes)
labels_colors(dend_RNAseq) <- colorCodes[groupCodes][order.dendrogram(dend_RNAseq)]
labels_colors(dend_protein) <- colorCodes[groupCodes][order.dendrogram(dend_protein)]


pdf("./outcome/clust_RNAseq_protein.pdf", height = 8, width = 6)
dendlist(dend_RNAseq, dend_protein) %>%
  untangle(method = "step2side") %>%
  tanglegram(
    common_subtrees_color_lines = FALSE, 
    highlight_distinct_edges  = FALSE, 
    highlight_branches_lwd=FALSE, 
    margin_inner = 8, margin_top = 2,
    margin_bottom = 3, margin_outer = 2,
    lwd=2, fast=TRUE
  )
dev.off()

## Differential expression analysis ------------------------------------------------
# RNAseq
library(tidyverse)
# There are some difference between the gene symbol form mapIds outcome and gene name in tx2gene file (generated from transcriptome reference)
tx2gene <- read.table("./data/RNAseq_data/tx2gene.txt", sep = "\t") %>%
  rename(transcriptID = V1) %>% rename(geneID = V2) %>% rename(geneName = V3)

# prepare files for pathway analysis in http://cpdb.molgen.mpg.de/ 
for (dose in c("Rif_The", "Rif_Tox")) {
  input_DEgenes <- read.table(paste0("./data/RNAseq_data/R-ODAF_test_RIF_DESeq2_RNA-Seq_",
                                     dose, "_vs_ConDMSO_FDR_0.01_DEG_table.txt")) %>% 
    rownames_to_column("geneID") %>% left_join(tx2gene[,c("geneID", "geneName")]) %>% 
    filter(padj < 0.01) %>% select(c("geneName")) %>%
    distinct() %>% 
    write.table(., paste0("./outcome/RNAseq_", dose, "_DEgeneNames.txt"),
                quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
}

# make pathway graph
RA_outcome <- list()
for(dose in c("Rif_The", "Rif_Tox")) {
  RA_outcome_temp <- read.csv(paste0("./outcome/RNAseq_", dose, "_ORA_results.tab"),
                              sep = "\t") %>% slice_head(n = 20) # select top pathways
  for (i in 1:nrow(RA_outcome_temp)) {
    RA_outcome_temp$genes[i] <- length(unlist(strsplit(RA_outcome_temp$members_input_overlap[i], ";")))
  }
  RA_outcome[[dose]] <- RA_outcome_temp
}

RNAseq_RAs <- RA_outcome$Rif_The[,c("pathway", "p.value", "genes")] %>% mutate(condition ="Rif_The") %>%
  full_join(RA_outcome$Rif_Tox[,c("pathway", "p.value", "genes")] %>% mutate(condition ="Rif_Tox"))
pdf("./outcome/RNAseq_RAs.pdf", height = 6, width = 8)
ggplot(data = RNAseq_RAs) + 
  geom_point(mapping = aes(x=-log(p.value,10), y=pathway, size = genes, colour = factor(condition))) +
  xlab("-log10(p.value)") + 
  xlim(0, max(plyr::round_any(-log(RNAseq_RAs$p.value,10), 10, f = ceiling)))
dev.off()

# get gene types
library(biomaRt)
listMarts()
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

RNAseq_DEgenetypes <- list()
for(dose in c("Rif_The", "Rif_Tox")) {
  RA_outcome_temp <- read.csv(paste0("./outcome/RNAseq_", dose, "_ORA_results.tab"), sep = "\t")
  gene_list_temp <- unique(gsub(" ", "", 
                              unlist(strsplit(RA_outcome_temp$members_input_overlap, ";")),
                              fixed = TRUE))
  annot_temp  <- getBM(attributes = c("hgnc_symbol", "gene_biotype", "description"),
                     filters = "hgnc_symbol",
                     values = read.table(paste0("./outcome/RNAseq_", dose, "_DEgeneNames.txt"), sep = "\t"),
                     mart = ensembl)
  RNAseq_DEgenetypes[[dose]] <- annot_temp %>% 
    mutate(inRAs = hgnc_symbol %in% gene_list_temp)
}
# check the information with Rif_The, can change the code for the Rif_Tox
RNAseq_DEgenetypes$Rif_The %>%
  group_by(gene_biotype) %>% summarise(n=n()) %>% mutate(prop=n/sum(n))
RNAseq_DEgenetypes$Rif_The %>% filter(gene_biotype == "lncRNA")

RNAseq_DEgenetypes$Rif_The %>%
  group_by(inRAs) %>% summarise(n=n()) %>% mutate(prop=n/sum(n))

# consider code to do (stop here)
which(duplicated(annot$hgnc_symbol)) # need to check duplicate
annot <- annot[!duplicated(annot$ensembl_gene_id),] # clean the duplicated gene Ids
annot %>% group_by(gene_biotype) %>% summarise(n = n()) %>% 
  ggplot(aes(reorder(gene_biotype, n), n)) + geom_bar(stat="identity") + xlab("type") + theme_bw() + coord_flip()

# Protein -------------------------------------------
load("./data/Rif_Protein_data.RData")
library(tidyverse)
View(ttest_results$The)
View(ttest_results$Tox)

protein_DEoverlap <- list()
for (dose in c("The", "Tox")) {
  protein_DEoverlap[[dose]] <- ttest_results$The %>% 
    dplyr::select(.,grep("qvalue", colnames(.))) %>% drop_na() %>%
    dplyr::filter_all(all_vars(. < 0.05))
}
# all proteins are in both Rif_The and Rif_Tox: 11 protein with "qvalue"
identical(rownames(protein_DEoverlap$The), rownames(protein_DEoverlap$Tox))
write.table(rownames(protein_DEoverlap$The),
            "./outcome/Protein_DEproteinIDs.txt",
            quote=FALSE, col.names=FALSE, row.names=FALSE, sep="\t")
# pathway analysis graph
RA_outcome <- read.csv("./outcome/Protein_Rif_qvalue_ORA_results.tab", 
                            sep = "\t") %>% slice_head(n = 20) # select the top pathways
for (i in 1:nrow(RA_outcome)) {
  RA_outcome$genes[i] <- length(unlist(strsplit(RA_outcome$members_input_overlap[i], ";")))
}

Protein_RAs <- RA_outcome[,c("pathway", "p.value", "genes")]
pdf("./outcome/Protein_RAs.pdf", height = 4.5, width = 7)
ggplot(data = Protein_RAs) + 
  geom_point(mapping = aes(x=-log(p.value,10), y=pathway, size = genes)) +
  xlab("-log10(p.value)") + 
  xlim(0, max(plyr::round_any(-log(Protein_RAs$p.value,10), 10, f = ceiling)))
dev.off()

## mape RNAseq to protein data -----------------------------------
protein_status <- ttest_results
names(protein_status) <- c("Rif_The", "Rif_Tox")

library(biomaRt)
listMarts()
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

RNAseq_DEgene_toProtein <- list()
for(dose in c("Rif_The", "Rif_Tox")) {
  RNAseq_temp  <- getBM(attributes = c("hgnc_symbol", "uniprotswissprot"),
                       filters = "hgnc_symbol",
                       values = read.table(paste0("./outcome/RNAseq_", dose,
                                                  "_DEgeneNames.txt"), sep = "\t"),
                       mart = ensembl)
  protein_dat <- protein_status[[dose]] %>% 
    dplyr::select(.,grep("qvalue", colnames(.))) %>% 
    filter(if_any(everything(), ~ !is.na(.)))
  
  RNAseq_DEgene_toProtein[[dose]] <- RNAseq_temp %>% 
    mutate(inProtein = uniprotswissprot %in% rownames(protein_dat)) %>%
    filter(inProtein ==TRUE) %>%
    write.table(., paste0("./outcome/RNAseq_DEgene_toProtein_", dose, "_geneNames.txt"),
              quote=FALSE, col.names=c("SYMBOL", "UNIPROT", "inProtein"), 
              row.names=FALSE, sep="\t")
}

DEgenes_inProtein <- cbind("Rif_The" = c(nrow(RNAseq_DEgene_toProtein$Rif_The),
                                         nrow(read.table("./outcome/RNAseq_Rif_The_DEgeneNames.txt"))),
                           "Rif_Tox" = c(nrow(RNAseq_DEgene_toProtein$Rif_Tox),
                                         nrow(read.table("./outcome/RNAseq_Rif_Tox_DEgeneNames.txt"))))
rownames(DEgenes_inProtein) <- c("DE genes detected at the protein level",
                                 "DE genes")
pdf("./outcome/DEgenes_inProtein.pdf", height = 4, width = 10)
barplot(DEgenes_inProtein, xlab = "Number of entities",
        col = c( "#F8766D", "#619CFF"), xlim = c(0, 6000),
        legend.text = rownames(DEgenes_inProtein),
        args.legend = list(x = "topright", inset = c(0, -0.3)),
        beside = FALSE,  horiz = TRUE)
dev.off()