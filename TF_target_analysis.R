# pipeline - link codes together (Need to move the folder later)
for (i in list.files("D:/Work/GitHub/RegmuNet/code")){
  source(paste0("D:/Work/GitHub/RegmuNet/code/", i))
}

sample_files <- list.files("./sample")
print(paste0("Number of sample files: ", length(sample_files), " files"))

metadat <- as.data.frame(matrix(NA, nrow = length(sample_files), ncol = 4))
colnames(metadat) <- c("file", "input_types",
                       "total_genes/proteins",
                       "annotation_status")
metadat$file <- sample_files 
metadat_v2 <- load_samples(metadat)
annotate(metadat_v2)
metadat_v2$input_types <- c(rep("Protein", 2), rep("RNAseq", 2))
TF_target_reference <- regula_database()
#TF_target_reference <- regula_database(tissue = "liver")

metadat_v3 <- get.exist_regula(metadat_v2, TF_target_reference)
metadat_v4 <- get.traceback_regula(metadat_v3, TF_target_reference)
metadat_v5 <- get.frac_traceback_regula(metadat_v4, TF_target_reference, cutoff = 0.05)
metadat_v6 <- across_TFs(metadat_v5) 

write.table(metadat_v6, "./TF_target_metadata.txt", quote=FALSE,
            col.names=TRUE, row.names=TRUE, sep="\t")
# sum up the TF-target analysis result -----------------------------------
get.merged_files <- function(files_list, merged_label) {
  library(tidyverse)
  input <- list()
  for (file_name in files_list) {
    input[[file_name]] <- read.table(file_name, header = TRUE) %>% 
      mutate({{file_name}} := TRUE)
  }
  outcome <- input %>% purrr::reduce(dplyr::full_join, by = merged_label) %>% 
    rename_with(~gsub("./sample/", "", .x)) %>%
    rename_with(~gsub("./output/", "", .x)) %>%
    rename_with(~gsub(".txt", "", .x))
  return(outcome)
}

setwd("D:/Work/GitHub/Rifampicin_omics_integration/")
files_list_1 <- c("./sample/annot_RNAseq_Rif_The_DEgeneNames.txt",
                "./sample/annot_RNAseq_Rif_Tox_DEgeneNames.txt",
                "./sample/annot_RNAseq_DEgene_toProtein_Rif_The_geneNames.txt",
                "./sample/annot_RNAseq_DEgene_toProtein_Rif_Tox_geneNames.txt",
                "./output/traceback_across_TFs.txt",
                "./output/existed_across_TFs.txt")

files_list_2 <- c("./output/traceback_TFs_RNAseq_Rif_The_DEgeneNames.txt",
                  "./output/traceback_TFs_RNAseq_Rif_Tox_DEgeneNames.txt",
                  "./output/traceback_TFs_RNAseq_DEgene_toProtein_Rif_The_geneNames.txt",
                  "./output/traceback_TFs_RNAseq_DEgene_toProtein_Rif_Tox_geneNames.txt")

outcome_TFs_info <- get.merged_files(files_list_1, merged_label = "SYMBOL") %>%
  full_join(get.merged_files(files_list_2, merged_label = "TF_traceback"),
            by = c("SYMBOL" = "TF_traceback")) %>% 
  dplyr::select(-starts_with("num_regula_targets")) %>%
  rename_with(~gsub("annot_", "", .x)) %>% 
  rename_with(~gsub("Names", "", .x)) %>%
  rename_with(~gsub("_gene", "", .x)) %>% 
  dplyr::rename(source = SYMBOL) %>% # change the name for the Cystoscape
  write.table(., "./outcome/outcome_TFs_info.txt", quote=FALSE,
              col.names=TRUE, row.names=FALSE, sep="\t")

setwd("D:/Work/GitHub/Rifampicin_omics_integration/output")
outcome_TF_relations <- lapply(list.files(pattern = "_relations_"),
                               function(x)read.table(x, header = TRUE) %>% 
                                 dplyr::select( c("TF", "target"))) %>%
  purrr::reduce(dplyr::full_join) %>% distinct() %>% 
  dplyr::rename(source = TF) #%>% # change the name for the Cystoscape
  write.table(., "./outcome_TF_relations.txt", quote=FALSE,
              col.names=TRUE, row.names=FALSE, sep="\t")
dim(outcome_TF_relations)

## Check the TF-target relation in MeDIP data analysis outcome ---------------------
MeDIP_Rif_The <- read.table("./outcome/MeDIP_selectedGenes_Rif_The_20220407.csv", 
                            header = TRUE, sep = ",")
MeDIP_Rif_Tox <- read.table("./outcome/MeDIP_selectedGenes_Rif_Tox_20220407.csv", 
                            header = TRUE, sep = ",")
MeDIP_genes <- unique(c(MeDIP_Rif_The$annot.symbol, MeDIP_Rif_Tox$annot.symbol))

list_files <- c("existed_relations_RNAseq_Rif_The_DEgeneNames.txt",
                "existed_relations_RNAseq_Rif_Tox_DEgeneNames.txt",
                "existed_relations_RNAseq_DEgene_toProtein_Rif_The_geneNames.txt",
                "existed_relations_RNAseq_DEgene_toProtein_Rif_Tox_geneNames.txt",
                "traceback_relations_RNAseq_Rif_The_DEgeneNames.txt",
                "traceback_relations_RNAseq_Rif_Tox_DEgeneNames.txt",
                "traceback_relations_RNAseq_DEgene_toProtein_Rif_The_geneNames.txt",
                "traceback_relations_RNAseq_DEgene_toProtein_Rif_Tox_geneNames.txt")
for (file_name in list_files) {
  dat_temp <- read.table(paste0("./output/", file_name), header = TRUE)
  print(file_name)
  print(MeDIP_genes[which(MeDIP_genes %in% dat_temp$TF)])
  print(MeDIP_genes[which(MeDIP_genes %in% dat_temp$target)])
}

outcome_MeDIPgenes <- list()
for (file_name in list_files) {
  dat_temp <- read.table(paste0("./output/", file_name), header = TRUE)
  outcome_MeDIPgenes[[file_name]] <- dat_temp[which(dat_temp$TF %in% MeDIP_genes | dat_temp$target %in% MeDIP_genes),]
}

outcome_SMARCA4 <- list()
for (file_name in list_files) {
  dat_temp <- read.table(paste0("./output/", file_name), header = TRUE)
  outcome_SMARCA4[[file_name]] <- dat_temp[which(dat_temp$TF == "SMARCA4" | dat_temp$target == "SMARCA4"),]
}

SMARCA4 <- outcome_MeDIPgenes$existed_relations_RNAseq_Rif_The_DEgeneNames.txt %>%
  filter(target == "SMARCA4")
MAN1B1 <- outcome_MeDIPgenes$existed_relations_RNAseq_Rif_The_DEgeneNames.txt %>%
  filter(target == "MAN1B1")
POLE <- outcome_MeDIPgenes$existed_relations_RNAseq_Rif_The_DEgeneNames.txt %>%
  filter(target == "POLE")
SMARCA4 <- outcome_MeDIPgenes$existed_relations_RNAseq_Rif_The_DEgeneNames.txt %>%
  filter(target == "SMARCA4")
a1 <- intersect(MAN1B1$TF, SMARCA4$TF)
a2 <- intersect(intersect(MAN1B1$TF, SMARCA4$TF), POLE$TF)
a1[-which(a1 %in% a2)]
MAN1B1$TF[-which(MAN1B1$TF %in% a1)]
SMARCA4$TF[-which(SMARCA4$TF %in% a1)]
POLE$TF[-which(POLE$TF %in% a2)]

a<- intersect(read.table("./outcome/RNAseq_Rif_The_DEgeneNames.txt"),
                 read.table("./outcome/RNAseq_Rif_Tox_DEgeneNames.txt"))
## select information for the network ------------------
frac_TFs <- top_frac(read.table("./output/existed_relations_RNAseq_Rif_The_DEgeneNames.txt",
                                header=TRUE) %>% dplyr::select(TF, target) %>%
                       dplyr::add_count(TF, sort = TRUE) %>% dplyr::select(TF, n) %>% 
                       distinct() %>% dplyr::rename("num_regula_targets" = "n"),
                     0.05)
a <- read.table("./output/existed_relations_RNAseq_Rif_The_DEgeneNames.txt",
           header=TRUE) %>% dplyr::select(TF, target) %>%
  dplyr::add_count(TF, sort = TRUE) %>% dplyr::select(TF, n) %>% 
  distinct() %>% dplyr::rename("num_regula_targets" = "n")
a_2 <- top_frac(a, 0.05)
b <- read.table("./output/existed_relations_RNAseq_Rif_Tox_DEgeneNames.txt",
                header=TRUE) %>% dplyr::select(TF, target) %>%
  dplyr::add_count(TF, sort = TRUE) %>% dplyr::select(TF, n) %>% 
  distinct() %>% dplyr::rename("num_regula_targets" = "n")
b_2 <- top_frac(b, 0.05)

selected_files <- c("./output/existed_relations_RNAseq_DEgene_toProtein_Rif_The_geneNames.txt")
selected_TF_relations <- lapply(list.files(pattern = "_relations_"),
                               function(x)read.table(x, header = TRUE) %>% 
                                 dplyr::select( c("TF", "target"))) %>%
  purrr::reduce(dplyr::full_join) %>% distinct() %>% 
  dplyr::rename(source = TF) #%>% # change the name for the Cystoscape
write.table(., "./outcome_TF_relations.txt", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")
dim(outcome_TF_relations)

# network ----------------------------------
library(igraph) # better for large input --> net to check the network analysis
library(tidyverse)

setwd("D:/Work/GitHub/Rifampicin_omics_integration")

nodes <- as_tibble(read.table("./outcome/outcome_TFs_info.txt", header = TRUE)) %>% 
  tibble::rowid_to_column("id")

edges <- as_tibble(read.table("./output/outcome_TF_relations.txt", header = TRUE)) %>% 
  left_join(nodes[, c("id", "source")], by = c("source")) %>%
  relocate(id) %>% dplyr::rename(from = id) %>% 
  left_join(nodes[, c("id", "source")], by = c("target" = "source")) %>%
  relocate(id, .after = target) %>% dplyr::rename(to = id) %>% 
  drop_na() %>% distinct()

selected_nodes <- nodes[which(nodes$source %in% edges$source & nodes$source %in% edges$target),]
selected_edges <- as_tibble(read.table("./output/outcome_TF_relations.txt", header = TRUE)) %>% 
  left_join(selected_nodes[, c("id", "source")], by = c("source")) %>%
  relocate(id) %>% dplyr::rename(from = id) %>% 
  left_join(selected_nodes[, c("id", "source")], by = c("target" = "source")) %>%
  relocate(id, .after = target) %>% dplyr::rename(to = id) %>% 
  drop_na() %>% distinct()

write.table(selected_nodes[,-1], "./network_selected_TFs.txt", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")
write.table(selected_edges[,c(2,3)], "./network_selected_TF_relations.txt", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")

write.table(full_join(selected_edges[,c(2,3)], selected_nodes[,-1]),
            "./network_selection.txt", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")
# focus on DDX5 - existed in both RNAseq and protein:
selected_net <-read.table("./outcome/network_selection.txt", header = TRUE)
selected_net_DDX5 <- selected_net[which(selected_net$source == "DDX5" | selected_net$target == "DDX5"),]
write.table(selected_net_DDX5,
            "./outcome/network_selection_forDDX5.txt", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")
source <- unique(c(selected_net_DDX5$source, selected_net_DDX5$target))

selected_TF_info <- as.data.frame(source) %>%
  left_join(read.table("./outcome/outcome_TFs_info.txt", header = TRUE)) %>%
  mutate(TF = ifelse(source %in% TF_target_reference$TF, TRUE, FALSE))


selected_TF_info <- as.data.frame(source) %>%
  left_join(read.table("./outcome/outcome_TFs_info.txt", header = TRUE)) %>%
  mutate(TF = ifelse(source %in% TF_target_reference$TF, TRUE, FALSE)) %>% 
  full_join(selected_net_DDX5) %>% relocate(target, .after = source) %>%
  mutate(DEgene = ifelse(RNAseq_Rif_The_DEgene==TRUE & RNAseq_Rif_Tox_DEgene == TRUE,
                         "both_conditions", "Rif_The/Tox")) %>%
  mutate(DEgene_inProtein = ifelse(RNAseq_DEgene_toProtein_Rif_The==TRUE & RNAseq_DEgene_toProtein_Rif_Tox== TRUE,
                         "both_conditions", "Rif_The/Tox")) %>%
  select(source, target, DEgene, DEgene_inProtein, TF)
write.table(selected_TF_info ,
            "./outcome/network_selection_forDDX5_info.txt", quote=FALSE,
            col.names=TRUE, row.names=FALSE, sep="\t")


#network_info <- list("nodes" = nodes, "edges" = edges)
network_info <- list("nodes" = selected_nodes, "edges" = selected_edges)

edges.igraph <- network_info$edges %>% dplyr::select(from, to)
nodes.igraph <-  network_info$nodes %>% rename(source = "label")
dat_igraph <- graph_from_data_frame(d = edges.igraph, vertices= nodes.igraph, 
                                    directed = TRUE)
print(dat_igraph, e=TRUE, v=TRUE)
pdf("test.pdf")
plot(dat_igraph)
dev.off()
## check number of gene in the TF-target
a <- read.table("./sample/annot_RNAseq_Rif_The_DEgeneNames.txt", header = TRUE)
a <- read.table("./sample/annot_RNAseq_Rif_Tox_DEgeneNames.txt", header = TRUE)
d_3 <- a[which(a$SYMBOL %in% TF_target_reference$TF | a$SYMBOL %in% TF_target_reference$target),]

# additional analysis with TF-target network
The <- read.table("./output/existed_relations_RNAseq_Rif_The_DEgeneNames.txt", header = TRUE)
The <- The[-c(66710:66713),]
Tox <- read.table("./output/existed_relations_RNAseq_Rif_Tox_DEgeneNames.txt", header = TRUE)

identical_TFs <- intersect(unique(The$TF), unique(Tox$TF))
identical_targets <- intersect(unique(The$target), unique(Tox$target))
identical_relations <- Tox[which(Tox$TF %in% identical_TFs & Tox$target %in% identical_targets),]

The_v2 <- The[which(The$TF %in% identical_TFs),]
Tox_v2 <- Tox[which(Tox$TF %in% identical_TFs),]
identical_targets_v2 <- intersect(unique(The_v2$target), unique(Tox_v2$target))

The_v3 <- The[-which(The$TF %in% identical_TFs),]
Tox_v3 <- Tox[-which(Tox$TF %in% identical_TFs),]
identical_targets_v3 <- intersect(unique(The_v3$target), unique(Tox_v3$target))

# test with network analysis
library(igraph)
get_paths_by_length <- function(g, len) {
  sp <- shortest.paths(g)
  sp[lower.tri(sp,TRUE)] <- NA
  wp <- which(sp==len, arr.ind=TRUE)
  mapply(function(a,b) get.shortest.paths(g, a, b)$vpath, wp[,1], wp[,2])
}

g2 <- graph_from_edgelist( as.matrix(identical_relations[,c(1,2)]),directed = FALSE)
g2 <- graph_from_edgelist( as.matrix(Tox_v3[,c(1,2)]),directed = FALSE)
g2 <- graph_from_edgelist( as.matrix(The_v3[,c(1,2)]),directed = FALSE)
g2 <- graph_from_edgelist( as.matrix(Tox[,c(1,2)]),directed = FALSE)
g2 <- graph_from_edgelist( as.matrix(The[,c(1,2)]),directed = FALSE)
diameter(g2)
sp<- shortest.paths(g2)
table(sp[upper.tri(sp)])
a<-get_paths_by_length(g2,4)

# track from proteomics
library(tidyverse)
The_TFs_pro <- top_frac(read.table("./output/traceback_TFs_RNAseq_DEgene_toProtein_Rif_The_geneNames.txt",
                         header = TRUE), 0.05)
The_v4 <- The[which(The$TF %in% c("DDX5", The_TFs_pro$TF_traceback)),]

Tox_TFs_pro <- top_frac(read.table("./output/traceback_TFs_RNAseq_DEgene_toProtein_Rif_Tox_geneNames.txt",
                                   header = TRUE), 0.05)
Tox_v4 <- Tox[which(Tox$TF %in% c("DDX5", Tox_TFs_pro$TF_traceback)),]

g2 <- graph_from_edgelist( as.matrix(Tox_v4[,c(1,2)]),directed = FALSE)
g2 <- graph_from_edgelist( as.matrix(The_v4[,c(1,2)]),directed = FALSE)
diameter(g2)
sp<- shortest.paths(g2)
table(sp[upper.tri(sp)])
a<-get_paths_by_length(g2,4)

## make graph
xy <- layout.fruchterman.reingold(g2)
seeds <- identical_TFs
# Colors
d <- distances(g2, to=V(g2)[ name %in% seeds], mode="out")
x <- apply(d, 1, function(x) paste(sort(colnames(d)[x == min(x)]), collapse="+"))
pal <- RColorBrewer::brewer.pal(length(unique(x)), "Set1")
cols <- pal[ match(x, sort(unique(x)))]

# Plot
plot(g2, 
     vertex.size=ifelse( V(g2)$name %in% seeds, 10, 5),
     vertex.color=cols,
     layout=xy, 
     vertex.label.color="black",
     # vertex.frame.color = par("bg"),
     edge.arrow.size=.5, 
     vertex.label=NA,
     margin=0
)
legend("bottomright",
       title="Dependencies",
       sort(unique(x)),
       pch=21,
       pt.bg=pal,
       col="black",
       bty="n"
)

# only select TFs involved
Tox <- read.table("./output/existed_relations_RNAseq_Rif_Tox_DEgeneNames.txt", header = TRUE)
f <- Tox[which(Tox$target %in% Tox$TF),c(1:2)]
g2 <- graph_from_edgelist( as.matrix(f[,c(1,2)]),directed = FALSE)
diameter(g2)
condition <-  c("normal", "normal", "check")
test_dat <- cbind(Tox[,c(1:2)],condition)
a<-graph.data.frame(test_dat)
diameter(a)
g2<-a
xy <- layout.fruchterman.reingold(g2)
seeds <- identical_TFs
# Colors
d <- distances(g2, to=V(g2)[ name %in% identical_TFs], mode="out")
x <- apply(d, 1, function(x) paste(sort(colnames(d)[x == min(x)]), collapse="+"))
pal <- RColorBrewer::brewer.pal(length(unique(test_dat$condition)), "Set1")
cols <- pal[ match(x, sort(unique(x)))]

vertex_attr(g2, "condition")
# Plot
plot(g2, 
     vertex.size=ifelse( V(g2)$name %in% c("check"), 10, 5),
     vertex.color=cols,
     layout=xy, 
     vertex.label.color="black",
     # vertex.frame.color = par("bg"),
     edge.arrow.size=.5, 
     vertex.label=NA,
     margin=0
)
legend("bottomright",
       title="Dependencies",
       sort(unique(x)),
       pch=21,
       pt.bg=pal,
       col="black",
       bty="n"
)

