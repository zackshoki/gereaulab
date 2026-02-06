library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(cowplot)

library(ExperimentHub)
eh <- ExperimentHub()
query(eh, "Kang")

(sce <- eh[["EH2259"]])

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
dim(sce)

# calculate per-cell quality control (QC) metrics
library(scater)
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]
dim(sce)

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)

# compute sum-factors & normalize
sce <- computeLibraryFactors(sce)
sce <- logNormCounts(sce)

sce$id <- paste0(sce$stim, sce$ind)
(sce <- prepSCE(sce, 
                kid = "cell", # subpopulation assignments
                gid = "stim",  # group IDs (ctrl/stim)
                sid = "id",   # sample IDs (ctrl/stim.1234)
                drop = TRUE))  # drop all other colData columns

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

# nb. of cells per cluster-sample

t(table(sce$cluster_id, sce$sample_id))
# compute UMAP using 1st 20 PCs
sce <- runUMAP(sce, pca = 20)
# wrapper to prettify reduced dimension plots
.plot_dr <- function(sce, dr, col)
  plotReducedDim(sce, dimred = dr, colour_by = col) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  theme_minimal() + theme(aspect.ratio = 1)

# downsample to max. 100 cells per cluster
cs_by_k <- split(colnames(sce), sce$cluster_id)
cs100 <- unlist(sapply(cs_by_k, function(u) 
  sample(u, min(length(u), 100))))

# plot t-SNE & UMAP colored by cluster & group ID
for (dr in c("TSNE", "UMAP"))
  for (col in c("cluster_id", "group_id"))
    plot5 <- .plot_dr(sce[, cs100], dr, col)
    plot5

    pb <- aggregateData(sce,
                        assay = "counts", fun = "sum",
                        by = c("cluster_id", "sample_id"))
    # one sheet per subpopulation
    assayNames(pb)
    t(head(assay(pb)))
    
    (pb_mds <- pbMDS(pb))
    # use very distinctive shaping of groups & change cluster colors
    pb_mds <- pb_mds + 
      scale_shape_manual(values = c(17, 4)) +
      scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2"))
    # change point size & alpha
    pb_mds$layers[[1]]$aes_params$size <- 5
    pb_mds$layers[[1]]$aes_params$alpha <- 0.6
    pb_mds
    
    # run DS analysis
    res <- pbDS(pb, verbose = FALSE)
    # access results table for 1st comparison
    tbl <- res$table[[1]]
    # one data.frame per cluster
    names(tbl)
    
    # view results for 1st cluster
    k1 <- tbl[[1]]
    head(format(k1[, -ncol(k1)], digits = 2))
    
    # construct design & contrast matrix
    ei <- metadata(sce)$experiment_info
    mm <- model.matrix(~ 0 + ei$group_id)
    dimnames(mm) <- list(ei$sample_id, levels(ei$group_id))
    contrast <- makeContrasts("stim-ctrl", levels = mm)
    
    # run DS analysis
    pbDS(pb, design = mm, contrast = contrast)
    
    
    # filter FDR < 5%, abs(logFC) > 1 & sort by adj. p-value
    tbl_fil <- lapply(tbl, function(u) {
      u <- dplyr::filter(u, p_adj.loc < 0.05, abs(logFC) > 1)
      dplyr::arrange(u, p_adj.loc)
    })
    
    # nb. of DS genes & % of total by cluster
    n_de <- vapply(tbl_fil, nrow, numeric(1))
    p_de <- format(n_de / nrow(sce) * 100, digits = 3)
    data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)
    frq <- calcExprFreqs(sce, assay = "counts", th = 0)
    # one sheet per cluster
    assayNames(frq)
    
    # expression frequencies in each
    # sample & group; 1st cluster
    t(head(assay(frq), 5))
    
    gids <- levels(sce$group_id)
    frq10 <- vapply(as.list(assays(frq)), 
                    function(u) apply(u[, gids] > 0.1, 1, any), 
                    logical(nrow(sce)))
    t(head(frq10))
    
    tbl_fil2 <- lapply(kids, function(k)
      dplyr::filter(tbl_fil[[k]], 
                    gene %in% names(which(frq10[, k]))))
    
    # nb. of DS genes & % of total by cluster
    n_de <- vapply(tbl_fil2, nrow, numeric(1))
    p_de <- format(n_de / nrow(sce) * 100, digits = 3)
    data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)
    
    # tidy format; attach pre-computed expression frequencies
    resDS(sce, res, bind = "row", frq = frq)
    
    # big-table (wide) format; attach CPMs
    resDS(sce, res, bind = "col", cpm = TRUE)
    
    library(UpSetR)
    de_gs_by_k <- map(tbl_fil, "gene")
    upset(fromList(de_gs_by_k))
    
    # pull top-8 DS genes across all clusters
    top8 <- bind_rows(tbl_fil) %>% 
      slice_min(p_adj.loc, n = 8, 
                with_ties = FALSE) %>% 
      pull("gene")
    
    # for ea. gene in 'top8', plot t-SNE colored by its expression 
    ps <- lapply(top8, function(g)
      .plot_dr(sce[, cs100], "TSNE", g) + 
        ggtitle(g) + theme(legend.position = "none"))
    
    # arrange plots
    plot_grid(plotlist = ps, ncol = 4, align = "vh")
    
    plotExpression(sce[, sce$cluster_id == "B cells"],
                   features = tbl_fil$`B cells`$gene[seq_len(6)],
                   x = "sample_id", colour_by = "group_id", ncol = 3) +
      guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # top-5 DS genes per cluster
    pbHeatmap(sce, res, top_n = 5)
    
    # top-20 DS genes for single cluster
    pbHeatmap(sce, res, k = "B cells")
    
    # single gene across all clusters
    pbHeatmap(sce, res, g = "ISG20")
