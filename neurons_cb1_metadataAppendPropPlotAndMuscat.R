# Zack Shoki 02/23/26 

rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(muscat)
library(scater)
library(lme4)      # dream depends on it
library(stats)
library(BiocParallel)
mem.maxVSize()
mem.maxVSize(32 * 1024^3)
register(MulticoreParam(workers = 4))
# tasks: 
# 1.) access the subsetted human neuron data
# 2.) add the metadata from "cann_metadata.csv" to that
# 3.) plot a proportion plot of the frequency of cannabis users to non users
# 4.) run muscat single-cell analysis


# 1.) 
setwd("/Users/zackshoki/Downloads/Shoki_Analysis") # adjust
neurons <- readRDS("neurons_CB1.rds")

cannabis_metadata <- read.csv("cann_metadata.csv")

# 2.)

# append a new slot to metadata (cannabis_use) of all zeros
neurons@meta.data$cannabis_use = ""

# iterate through each donor code that appears in the metadata (there might be an easier, logical indexing type way to do this)
for (donor_code in unique(neurons@meta.data$Donor)) {
  #   is the donor code contained in cannabis_metadata
  if (donor_code %in% cannabis_metadata$Donor) {
    use_cannabis <- cannabis_metadata[cannabis_metadata$Donor == donor_code, "cann_user"]
    #   if yes, 
    #     check for a "Y" in cannabis_metadata donor row
    if (use_cannabis == "Y") {
    #       if yes, 
    #         cannabis_use at the donor's row in neurons = "Y"
      neurons@meta.data[neurons@meta.data$Donor == donor_code, "cannabis_use"] <- "Y"

    } else {
      #   else, 
      #    cannabis_use at the donor's row in neurons = "N"
      neurons@meta.data[neurons@meta.data$Donor == donor_code, "cannabis_use"] <- "N"
    }
  } else {
    #   if no, 
    #      cannabis_use at the donor's row in neurons = "unknown"
    neurons@meta.data[neurons@meta.data$Donor == donor_code, "cannabis_use"] <- "unknown"
    
  }
  
}

# subset out cells with cannabis use of unknown
neurons <- subset(neurons, subset = neurons@meta.data$cannabis_use != "unknown")

# saveRDS(neurons, file = "neurons_cb1_withCannMetaData.RDS")

# 3.)

# create proportion dataframe
prop_df <- neurons@meta.data %>%
  group_by(cannabis_use) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

# create proportion barplot
plot1 <- ggplot(prop_df, aes(x = "", y = prop, fill = cannabis_use)) +
  geom_bar(stat = "identity", width = 0.55) +
  theme_minimal() +
  labs(fill = "", x = NULL, title = "Proportion of Cells from Cannabis Users to Cells from Non-Cannabis Users") +
  scale_fill_manual(values = c("N" = "lightblue", "Y" = "pink"), labels = c("N" = "Non-Cannabis User", "Y" = "Cannabis User")) +
  ylab("Proportion of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

# print plot
plot1

# create folder for plot
dir_name <- "neurons_cb1Plots"

if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

# save plot to folder
ggsave(
  filename = file.path(dir_name, paste0("CannabisUserPropPlot.pdf")),
  plot = plot1,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf
)

# 4.) 

# take out cells with unknown donors
neurons <- subset(neurons, subset = Donor != "unassigned" & Donor != "unk")

# create sce for cannabis users and non users
users <- as.SingleCellExperiment(subset(neurons, subset = cannabis_use == "Y"))
non_users <- as.SingleCellExperiment(subset(neurons, subset = cannabis_use == "N"))


# quality control
users <- users[rowSums(counts(users) > 0) > 0, ]

qc <- perCellQCMetrics(users)

ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
users <- users[, !ol]
users <- users[rowSums(counts(users) > 1) >= 10, ]

users <- computeLibraryFactors(users)
users <- logNormCounts(users)


non_users <- non_users[rowSums(counts(non_users) > 0) > 0, ]

qc <- perCellQCMetrics(non_users)

ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
non_users <- non_users[, !ol]
non_users <- non_users[rowSums(counts(non_users) > 1) >= 10, ]

non_users <- computeLibraryFactors(non_users)
non_users <- logNormCounts(non_users)

# prepare objects for muscat analysis

users <- prepSCE(
  users,
  kid = "Atlas_annotation",  # cluster assignments
  gid = "cannabis_use",    # experimental groups
  sid = "Donor",   # unique sample IDs
  drop = TRUE          # drop all other metadata columns
)

non_users <- prepSCE(
  non_users,
  kid = "Atlas_annotation",  # cluster assignments
  gid = "cannabis_use",    # experimental groups
  sid = "Donor",   # unique sample IDs
  drop = TRUE          # drop all other metadata columns
)

# combine users and non_users back to one
common_genes <- intersect(rownames(users), rownames(non_users))

users_sub   <- users[common_genes, ]
non_users_sub <- non_users[common_genes, ]

sce <- cbind(users_sub,  non_users_sub)

# define number of clusters, number of samples, names of clusters, and names of samples as variables
nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# compute UMAP using 1st 20 PCs
sce <- runUMAP(sce, pca = 20)
sce <- runTSNE(sce, pca = 20)

.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

# plot t-SNE & UMAP colored by cluster & group ID TODO: SAVE PLOTS
plots <- list()
idx <- 1
for (dr in c("UMAP", "TSNE")) {
  for (col in c("cluster_id", "group_id")) {
    plots[[idx]] <- .plot_dr(sce, dr, col) 
    idx <- idx + 1
  }
}
  

## check plots.. is there adequate representation from cannabis_users of all cell types?
    for (i in seq_along(plots)) {
      
      ggsave(
        filename = file.path(dir_name, paste0("plot_", i, ".pdf")),
        plot = plots[[i]],
        width = 20,
        height = 12,
        units = "in",
        device = cairo_pdf
      )
    }
    
mm <- mmDS(
  sce,
  method = "dream",
  BPPARAM = MulticoreParam(workers = 4),
  n_cells = 2,       # minimum cells per cluster per donor
  min_cells = 1      # minimum cells expressing a gene in cluster
)
