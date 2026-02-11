#020626
#Zack practice code
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("C:/Users/julietm/Box/Shoki_Analysis")
neurons <- readRDS("Seurat_shams_final.Rds")
genes <- c("Scn10a",
           "Pvalb",
           "Osmr",
           "Cnr1",
           "Nppb",
           "Sst",
           "Calca",
           "Mrgprd",
           "Ntrk3",
           "Trpm8",
           "Atf3",
           "S100a16",
           "Ntrk2",
           "Th",
           "Trpv1",
           "Oprk1",
           "Sstr2",
           "Adra2a",
           "Mrgprb4")

neurons
neurons@meta.data
unique(neurons@meta.data$Species)

human_neurons <- subset(
  neurons, 
  subset = Species == "Human"
)
mouse_neurons <- subset(
  neurons,
  subset = Species == "Mouse"
)
human_neurons
human_neurons@meta.data
unique(human_neurons@meta.data$Species)

mouse_neurons

# <- group_by(neurons) filter by species???
unique(neurons@meta.data$seurat_clusters)
unique(neurons@meta.data$final)

# finding cell markers
markers <- FindMarkers(neurons, ident.1 = "Trpm8")
# go into the dataset, find the cells that have identities matching the features arg, plot them
plot1 <- RidgePlot(neurons, features = genes[c(2, 9)])

# checking how neurons are grouped currently (by gene)
unique(Idents(neurons))
unique(genes)

levels(neurons@meta.data$final)

# notes:
# cells are clustered by Seurat algorithm into groups based on similarity, numbered (theres 32 of them here)
# clusters are grouped into groups (populations) based on a specific marker gene that is heavily expressed, naming scheme varies, but here, the marker gene = the name of the cell population (theres like 18 of em here)

# TODO: 
# filter the data set for each species (or sex, age, etc.)
# plot specific cell populations for different conditions
# group the clusters into populations with seurat's FindAllMarkers, research each gene that is highly expressed in a given population and different neuron types, then try to inductively label what each population type is

#020926

# lets compare genes of human females to human males

human_males <- subset(
  human_neurons,
  subset = Sex == "M"
)

human_females <- subset(
  human_neurons,
  subset = Sex == "F"
)

plot2 <- RidgePlot(human_neurons, features = genes[1:5], group.by = "Sex")
plota <- RidgePlot(human_females, features = genes[1:5])
plotb <- RidgePlot(human_males, features = genes[1:5])
plot3 <- plota+plotb
plot4 <- DimPlot(neurons, group.by = "Species")
plot5 <- FeaturePlot(neurons, features = "Pvalb")
plot6 <- DotPlot(neurons, features = genes)

avg <- AverageExpression(neurons, features = genes, group.by = "final")$RNA
peak_cluster <- apply(avg, 1, which.max)
genes_ordered <- names(sort(peak_cluster))

plot7 <- DotPlot(
  neurons,
  features = genes_ordered[1:5],
  group.by = "final"
)

# Now lets erase the clustering that was precreated and try to create new clusters 
# and see what each gene is actually doing by comparing expression 
# across different conditions

neurons.fresh <- neurons

neurons.fresh@meta.data$seurat_clusters <= NULL



neurons.fresh <- FindVariableFeatures(neurons.fresh, selection.method = "vst")

neurons.fresh <- ScaleData(neurons.fresh)

plot8 <- ElbowPlot(neurons.fresh)

neurons.fresh <- RunPCA(neurons.fresh, npcs = 20)

neurons.fresh <- FindNeighbors(neurons.fresh, dims = 1:20)

neurons.fresh <- FindClusters(neurons.fresh, resolution = 0.5)

plot9 <- RidgePlot(neurons.fresh, features = genes[1:5])
plot10 <- RidgePlot(neurons.fresh, features = genes[6:10])
plot11 <- RidgePlot(neurons.fresh, features = "Pvalb", group.by = "seurat_clusters")

# im thinking cluster 13 are the "Pvalb" cells, so proprioceptive function??  
plot12 <- DimPlot(neurons.fresh, group.by = "seurat_clusters")
# Pvalb seems may be closely associated with cluster 2
plot13 <- FeaturePlot(neurons.fresh, features = genes)
plot14 <- DotPlot(neurons.fresh, features = genes, group.by = "seurat_clusters")

# Pvalb seems to track as cluster 13, ntrk2 might be cluster 7(encodes a protein), osmr is likely cluster 4 (induces itch?)
# it seems that osmr expression may vary from species to species so we might check that soon

#fresh.markers <- FindAllMarkers(neurons.fresh, only.pos = TRUE)
#fresh.markers <- fresh.markers %>%
#  group_by(cluster) %>%
#  dplyr::filter(avg_log2FC > 1)

plots <- list(
  plot1,
  plot2,
  plot3,
  plot4,
  plot5, 
  plot6, 
  plot7, 
  plot8, 
  plot9, 
  plot10,
  plot11,
  plot12, 
  plot13,
  plot14
)
for (i in seq_along(plots)) {
  png(filename = paste0("plot_", i, ".png"), width = 800, height = 600)
  print(plots[[i]])     
  dev.off()      
}