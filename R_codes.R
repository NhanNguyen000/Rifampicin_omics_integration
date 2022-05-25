# preparing toy data -----------------------------------
# packages & parameters 
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
p_cutoff <- 0.05

# import protein data 
load("./data/Rif_Protein_data.RData")
Protein_norm_dat <- dat_sample
DEPs <- list()
for (dose in names(ttest_results)) {
  for (time in unique(substring(colnames(ttest_results[[dose]]), 9,11))) {
    dat_temp <- ttest_results[[dose]] %>% 
      dplyr::select(contains(time)) %>% rename_with(~substring(.x,13)) %>% 
      filter(pvalue < p_cutoff) %>% filter(qvalue < p_cutoff) %>%
      rownames_to_column(var="UNIPROT")
    dat_temp$SYMBOL <- mapIds(org.Hs.eg.db, keys = dat_temp$UNIPROT,
                              column = "SYMBOL", keytype = "UNIPROT", multiVals = "first") #need to check if rerutn the correct row
    DEPs[[paste0(dose,"_", time)]] <- dat_temp
  }
}
rm(dat_sample, ttest_results, dose, time)

# import RNAseq data (outcome from basic DESeq2 analysis)
library(DESeq2)
load("./data/Rif_RNAseq_DEseq2_allsamples.RData")
RNAseq_norm_dat <- assay(norm_data)
DEGs_allsamples <- as_tibble(DESeq_result) %>% filter(pvalue < p_cutoff) %>% filter(padj < p_cutoff)
DEGs_allsamples$ENSEMBL <- mapIds(org.Hs.eg.db, keys = DEGs_allsamples$SYMBOL,
                                  column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
RNAseq_norm_DEGs_allsamples <- as.data.frame(RNAseq_norm_dat) %>% 
  mutate(ENSEMBL=sub("\\..*", "", row.names(RNAseq_norm_dat))) %>% filter(ENSEMBL %in% DEGs_allsamples$ENSEMBL) %>%
  dplyr::select(-ENSEMBL)
rm(DESeq_outcome, DESeq_result, norm_data)

load("./data/Rif_RNAseq_DEseq2_perConditions.RData")
p_cutoff <- 0.01
names(outcome) <- substring(names(outcome),5,11)
DEGs <- list()
for (condition in names(outcome)) {
  DEGs[[condition]] <- as_tibble(outcome[[condition]]$DESeq_result) %>% filter(pvalue < p_cutoff) %>% filter(padj < p_cutoff)
}
rm(outcome)

# import RNAseq data (outcome from DESeq2 analysis + R-ODAF)
DEGs_roadf <- list()
folder <- "./data/HeCaTos_Hepatic/The_DMSO_ctrl/deseq2_output/" 
for (condition in c("The_002", "The_008", "The_024", "The_072")) {
  dat_temp <- read.csv(paste0(folder, "R-ODAF_", condition,
                              "_DESeq2_RIF_vs_DMSO_FDR_0.01_DEG_table.txt"),
                       sep = "\t")
  for(type in c("SYMBOL", "ENTREZID", "GENENAME")) {
    dat_temp[[type]] <- mapIds(org.Hs.eg.db,
                               keys = sub("\\..*", "", row.names(dat_temp)),
                               column = type, keytype = "ENSEMBL", multiVals = "first")
  }
  DEGs_roadf[[condition]] <- dat_temp
}
# import MeDIP data (waiting)

# overview about the improt data --------------------------
# plot the number of DEGs, DEPs in one chart
Dose <- substring(names(DEGs), 1,3)
Time <- substring(names(DEGs), 5,7)

DEGs_num <- c()
DEPs_num <- c()
for(i in names(DEGs)) {
  DEGs_num <- c(DEGs_num, nrow(DEGs[[i]]))
  DEPs_num <- c(DEPs_num, nrow(DEPs[[i]]))
}

dat_plot <- as.data.frame(cbind(Dose, Time, as.numeric(DEGs_num), as.numeric(DEPs_num)))

ggplot(dat_plot, aes(x=Time, y=DEGs_num, colour = Dose)) + 
  geom_point(aes(shape=Dose), size = 5) + geom_line(aes(group = Dose)) + 
  scale_y_continuous(limits=c(0, 6000)) +
  theme(text = element_text(size=15))

ggplot(dat_plot, aes(x=Time, y=DEPs_num, colour = Dose)) + 
  geom_point(aes(shape=Dose), size = 5) + geom_line(aes(group = Dose)) + 
  scale_y_continuous(limits=c(0, 500)) +
  theme(text = element_text(size=15))

# Make PCA for each omics 
library(factoextra)
library("viridis")   
dat <- t(RNAseq_norm_dat)
#dat <- t(RNAseq_norm_DEGs_allsamples) # not able to separate the Rif samples vs. control
res.pca <- prcomp(dat)
group <- as.factor(substring(rownames(dat), 1, 7))
fviz_pca_ind(res.pca, label = substring(rownames(dat), 9,11),
             addEllipses = TRUE, habillage=group)

# Should make the PCA nicer with tidyver, and also make PCA for the protein data

# Making gene regulatory network -------------------------------------------
geneTF_info <-read.csv("./data/hTFtarget_information.txt", sep = "\t")

# toy data 1 = option 1: make one big data table - select with rank >= 3
identical(names(DEGs), names(DEPs)) # TRUE
outcome <- list()
for (condition in names(DEGs)) {
  dat_gene <- as_tibble(DEGs[[condition]]) %>% 
    dplyr::select(gen_baseMean = baseMean, gen_log2FC = log2FoldChange,
                  gen_pvalue = pvalue, gen_padj = padj,
                  SYMBOL, ENTREZID, GENENAME)
  dat_pro <- DEPs[[condition]] %>% rownames_to_column(var="UNIPROT") %>%
    left_join(Protein_annot) %>% 
    dplyr::rename(pro_pvalue = pvalue, pro_qvalue = qvalue,
                  pro_effectSize = effectSize, pro_directedEffect = directedEffect)
  mega_gen_pro <- dat_gene %>% full_join(dat_pro, by=c("SYMBOL", "ENTREZID", "GENENAME")) %>%
    mutate(isTF = ifelse(SYMBOL %in% geneTF_info$TF, TRUE, FALSE)) %>%
    mutate(isTarget = ifelse(SYMBOL %in% geneTF_info$target, TRUE, FALSE)) %>%
    mutate(isDEG = ifelse(gen_padj < 0.05, TRUE, FALSE)) %>%
    mutate(isDEP = ifelse(pro_qvalue < 0.05, TRUE, FALSE)) %>%
    rowwise() %>% mutate(rank = sum(isTF, isTarget, isDEG, isDEP, na.rm = TRUE))
  outcome_temp <- mega_gen_pro %>% 
    dplyr::select(SYMBOL, ENTREZID, GENENAME, UNIPROT,
                  isTF, isTarget, isDEG, isDEP, rank) %>% filter(rank >= 3)
  outcome[[condition]] <- outcome_temp
  rm(dat_gene, dat_pro, mega_gen_pro, outcome_temp)
}
# toy data 2 = option 2: build a network form the begining
DEGs <- DEGs_roadf
outcome <- list()
for (condition in names(DEGs)) {
  dat_temp <- geneTF_info %>% filter(TF %in% DEGs[[condition]]$SYMBOL | TF %in% DEPs[[condition]]$SYMBOL) %>%
    filter(target %in% DEGs[[condition]]$SYMBOL | target %in% DEPs[[condition]]$SYMBOL) %>%
    mutate(TF_isDEG = ifelse(TF %in% DEGs[[condition]]$SYMBOL, TRUE, FALSE)) %>%
    mutate(TF_isDEP = ifelse(TF %in% DEPs[[condition]]$SYMBOL, TRUE, FALSE)) %>%
    mutate(target_isDEG = ifelse(target %in% DEGs[[condition]]$SYMBOL, TRUE, FALSE)) %>%
    mutate(target_isDEP = ifelse(target %in% DEPs[[condition]]$SYMBOL, TRUE, FALSE)) %>%
    rowwise() %>% mutate(weighted = sum(TF_isDEG, TF_isDEP, target_isDEG, target_isDEP, na.rm = TRUE))
  names(dat_temp)[names(dat_temp) == "weighted"]<- condition
  outcome[[condition]] <- dat_temp
}

test <- outcome %>% purrr::reduce(full_join, by = c("TF", "target", "tissue")) %>%
  dplyr::select(-starts_with("TF_")) %>% dplyr::select(-starts_with("target_")) %>%
  mutate(weight = sum(c_across(contains("_")), na.rm=TRUE)) 
test <- test %>% dplyr::filter(weight > 6) %>% dplyr::select(TF, target, weight)
test <- test %>% filter(weight >26) %>% dplyr::select(TF, target, weight)
# --> prepare for the network:
nodes <- data.frame("label" = unique(c(test$TF, test$target))) %>% tibble::rowid_to_column("id") %>%
  mutate(type = ifelse(label %in% geneTF_info$TF, "TF", "target"))
edges <- test %>% left_join(nodes, by = c("TF" = "label")) %>% left_join(nodes, by = c("target" = "label")) %>%
  dplyr::select(id.x, id.y, weight)
names(edges) <- c("from", "to", "weight")

# toy data 3: 
get.TFs <- function(Symbol_list, TF_target_dat) {
  TFs <- Symbol_list[Symbol_list %in% TF_target_dat$TF]
  return(TFs)
}
get.targets <- function(Symbol_list, TF_target_dat) {
  targets <- Symbol_list[Symbol_list %in% TF_target_dat$target]
  return(targets)
}

get.tabl <- function(SYMBOL, condition) {
  if (length(SYMBOL) >= 1) {
    outcome <- cbind(SYMBOL,"Boolean" = TRUE)
    colnames(outcome)[2] <- condition
  }
  if (length(SYMBOL) == 0) outcome <- NULL
  return(as.data.frame(outcome))
}

TF_target_dat <- geneTF_info
Selected_DEGs_TF <- list()
Selected_DEGs_target <- list()
for (condition in names(DEGs)) {
  TF_temp <- get.tabl(get.TFs(DEGs[[condition]]$SYMBOL, TF_target_dat),
                      condition)
  target_temp <- get.tabl(get.targets(DEGs[[condition]]$SYMBOL, TF_target_dat),
                          condition)

  if (length(TF_temp) >0) Selected_DEGs_TF[[condition]] <- TF_temp
  if (length(target_temp) >0) Selected_DEGs_target[[condition]] <- target_temp
}

Selected_DEPs_TF <- list()
Selected_DEPs_target <- list()
for (condition in names(DEPs)) {
  TF_temp <- get.tabl(get.TFs(DEPs[[condition]]$SYMBOL, TF_target_dat),
                      condition)
  target_temp <- get.tabl(get.targets(DEPs[[condition]]$SYMBOL, TF_target_dat),
                          condition)
  
  if (length(TF_temp) >0) Selected_DEPs_TF[[condition]] <- TF_temp
  if (length(target_temp) >0) Selected_DEPs_target[[condition]] <- target_temp
}

DEG_TFs <- Selected_DEGs_TF %>% purrr::reduce(full_join, by = "SYMBOL")
DEG_targets <- Selected_DEGs_target %>% purrr::reduce(full_join, by = "SYMBOL")
DEP_TFs <- Selected_DEPs_TF %>% purrr::reduce(full_join, by = "SYMBOL")
DEP_targets <- Selected_DEPs_target %>% purrr::reduce(full_join, by = "SYMBOL")

select_TFs <- unique(DEG_TFs$SYMBOL, DEP_TFs$SYMBOL)
select_targets <- unique(DEG_targets$SYMBOL, DEP_targets$SYMBOL)


select_TFs <- intersect(DEG_TFs$SYMBOL, DEP_TFs$SYMBOL)
select_targets <- intersect(DEG_targets$SYMBOL, DEP_targets$SYMBOL)

a<- geneTF_info[which(geneTF_info$TF==select_TFs),]
b<- geneTF_info[which(geneTF_info$target==select_TFs),]

b2 <- b[which(b$TF %in% DEG_TFs$SYMBOL),]
# --> prepare for the network:
nodes <- data.frame("label" = unique(c(select_TFs, select_targets ))) %>% 
  tibble::rowid_to_column("id") %>%
  mutate(type = ifelse(label %in% geneTF_info$TF, "TF", "target"))

edges <- nodes %>% left_join(TF_target_dat, by = c("label" = "TF")) %>%
  filter(target %in% select_targets) %>% dplyr::select(-c("tissue")) %>% 
  distinct() %>%
  left_join(nodes, by = c("target" = "label")) %>% dplyr::select(-starts_with("type"))
names(edges) <- c("from", "TF", "target", "to")
edges <- edges[,-c(2,3)]

# making interactive network
library(visNetwork) # give interactive + nice format network but for small input
nodes.vis <- nodes %>% mutate(group = type)
edges.vis <- edges %>% 
  mutate(arrows = "to") %>%
  mutate(color = "red")
#edges.vis$label <- paste0("all", c(1:14))
save(nodes.vis, edges.vis, file = "visNetwork.RData")
visnet <- visNetwork(nodes.vis, edges.vis) %>% 
  visGroups(groupname = "TF", color = "orange") %>%
  visGroups(groupname = "target", color = "lightblue") %>%
  visLegend() %>% visInteraction(navigationButtons = TRUE)
visOptions(visnet, 
           nodesIdSelection = TRUE, highlightNearest = TRUE,
           selectedBy = "group") %>%
  visSave(file = "network.html", background = "white")

# making static network - for Grow Science day 2021:
library(igraph) # better for large input --> net to check the network analysis
dat_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
plot(dat_igraph, edge.arrow.size=.4, edge.curved=.1)

deg <- degree(dat_igraph, mode="all")
hist(deg, main="Histogram of node degree")
deg.dist <- degree_distribution(dat_igraph, cumulative=T, mode="all")
plot(x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
     xlab="Degree", ylab="Cumulative Frequency")
mean_distance(dat_igraph, directed=F)
range(distances(dat_igraph))
ceb <- cluster_edge_betweenness(dat_igraph, directed = F) 
dendPlot(ceb, mode="hclust")
V(dat_igraph)$label.cex <- 0.6
pdf(file = "GrowScienceDay_2021.pdf", width = 15, height = 8)
plot(ceb, dat_igraph,
     edge.arrow.size=.4, edge.curved=.2) 
dev.off()
