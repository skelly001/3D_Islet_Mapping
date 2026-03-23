# ==============================================================================
# 3D Spatial Proteomics Analysis with Giotto
# ==============================================================================
# Script: 0a_1-3D_spatial_proteomics.R
# Description: Performs 3D spatial proteomics analysis using the Giotto
#              framework. Includes data preprocessing (filtering, normalization),
#              statistics and spatial visualization, HVG feature selection,
#              dimensionality reduction (PCA, UMAP), Leiden clustering,
#              differential expression, and Leiden cluster metagene
#              creation and visualization.
#
# Input:  - {data_directory}/expression.csv
#         - {data_directory}/spatial_locs.csv
# Output: - Filter combination plot
#         - 2D/3D spatial plots (nr_genes)
#         - PCA plots (unlabeled and Leiden-colored)
#         - Metadata heatmap (top DEGs per Leiden cluster)
#         - Spatial metagene plots (Leiden cluster metagenes)
#           (all plots saved to save_directory via Giotto instructions)
# ==============================================================================



library(tidyverse)
library(ggplot2)
library(reshape2)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggpubr)
remotes::install_github("RubD/Giotto")
library(Giotto)


#1_Environment setup
installGiottoEnvironment()


#2_Downloading dataset
data_directory <- "data/0-3D_spatial_proteomics"
save_directory <- "output/RD0-3D_spatial_proteomics"


#3_Creating a Giotto object
my_instructions <- createGiottoInstructions(save_plot = TRUE,
                                            show_plot = TRUE,
                                            return_plot = FALSE,
                                            save_dir = save_directory)

my_giotto_object<-createGiottoObject(raw_exprs=paste0(data_directory,
                                                      "/expression.csv"),
                                     spatial_locs = paste0(data_directory,
                                                           "/spatial_locs.csv"),
                                     instruction=my_instructions)


#4_Preprocessing
##4-1_Filtering
filterCombinations(gobject = my_giotto_object,
                   expression_thresholds = c(1E+4, 1E+5),
                   gene_det_in_min_cells = c(3, 5, 10),
                   min_det_genes_per_cell = c(1000, 2000, 3000))

##expression_thresholds	;all thresholds to consider a gene expressed
##gene_det_in_min_cells	;minimum number of cells that should express a gene to consider that gene further
##min_det_genes_per_cell	;minimum number of expressed genes per cell to consider that cell further


my_giotto_object <- filterGiotto(gobject = my_giotto_object,
                                 expression_threshold = 3,
                                 gene_det_in_min_cells = 3,
                                 min_det_genes_per_cell = 100)


#4-2_Normalization
my_giotto_object <- normalizeGiotto(gobject = my_giotto_object,
                                    norm_methods = "standard",
                                    scalefactor = 6000,
                                    scale_order = "first_genes")

#4-3_Statistics
my_giotto_object <- addStatistics(gobject = my_giotto_object)
# view gene and cell stats respectively
head(fDataDT(my_giotto_object))
head(pDataDT(my_giotto_object))

###number of genes
spatPlot2D(gobject = my_giotto_object,
           show_image = TRUE,
           point_alpha = 1,
           point_size = 5,
           cell_color = 'nr_genes',
           color_as_factor = F)

spatPlot3D(gobject = my_giotto_object,
           axis_scale = "cube",
           point_size = 5,
           cell_color = 'nr_genes')



#5_Clustereing and cell-type identification #HVG

##5-1_Feature selection
my_giotto_object <- calculateHVG(gobject = my_giotto_object,
                                 expression_values = "normalized",
                                 method = "cov_groups",
                                 nr_expression_groups = 20,
                                 zscore_threshold = 1.5)

##5-2_Dimensionality reduction
my_giotto_object <- runPCA(gobject = my_giotto_object,
                           expression_values = "normalized",
                           genes_to_use = "hvg")



plotPCA(gobject = my_giotto_object,
        point_size=5,
        show_legend=TRUE,
        show_plot=TRUE,
        legend_text=10)

my_giotto_object<-runUMAP(gobject = my_giotto_object,
                          dimensions_to_use = 1:10,
                          n_neighbors = 40,
                          n_components = 2,
                          min_dist = 0.01)



##5-3_Leiden cluster
my_giotto_object <- doLeidenCluster(gobject = my_giotto_object,
                                    name = "leiden_clus")


plotPCA(gobject = my_giotto_object,
        point_size=5,
        show_legend=TRUE,
        show_plot=TRUE,
        legend_text=10,
        cell_color='leiden_clus',
        color_as_factor = TRUE)




##5-5_DEGs, Heatmap
ST_scran_markers_subclusters = findMarkers_one_vs_all(gobject = my_giotto_object,
                                                      method = 'scran',
                                                      expression_values ='normalized',
                                                      cluster_column = 'leiden_clus')

ST_top3genes = ST_scran_markers_subclusters[, head(.SD, 3), by = 'cluster']$genes
ST_topNgenes = ST_scran_markers_subclusters[, head(.SD, 30), by = 'cluster']$genes

plotMetaDataHeatmap(gobject = my_giotto_object,
                    selected_genes = ST_topNgenes,
                    metadata_cols = c('leiden_clus'),
                    gradient_color = c("#440154", "#21908C", "#FDE725"),
                    custom_cluster_order = c("1", "3", "4", "2"))


##############################################################
#Create metagenes from cluster modules and visualize_Leiden###
###############################################################

top40_per_leiden=cluster_genes_DT[,head(.SD, 40), by = clus]
topN_per_leiden=ST_scran_markers_subclusters[,head(.SD,40), by='cluster']

Leiden_cluster_genes=topN_per_leiden$cluster; names(Leiden_cluster_genes)= topN_per_leiden$genes

my_giotto_object= createMetagenes(my_giotto_object,
                                  gene_clusters = Leiden_cluster_genes,
                                  name = 'Leiden_cluster_metagene')


spatCellPlot(my_giotto_object,
             spat_enr_names='Leiden_cluster_metagene',
             cell_annotation_values=as.character(c(1:8)),
             point_size=5,
             cow_n_col=2)
