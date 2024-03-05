#from limoR
setwd("/home/guanpengpeng/jupyterlab/Brain_remove_double_mt2/try_sample3/num14_cell_type/miloR")

###through sceasy, the h5ad format can be converted to rds
library(sceasy)
Sys.setenv(RETICULATE_PYTHON = "/home/guanpengpeng/software/Anaconda/install/bin/python")
#h5ad_file = "../Neuro_umap_bbknn.OPCs.disease_diff.h5ad"
h5ad_file = "num2.h5ad"
sceasy::convertFormat(h5ad_file, from="anndata", to="seurat",outFile='num2.rds')

###run milo
#miloR github https://github.com/MarioniLab/miloR
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

###test the miloR using num2.rds
GSE138852_data <- readRDS("./GSE138852_test/GSE138852_filter.rds")
GSE138852_sce <- as.SingleCellExperiment(GSE138852_data)

GSE138852_sce_milo <- Milo(GSE138852_sce)
GSE138852_sce_milo <- buildGraph(GSE138852_sce_milo, k = 100, d = 50)
#num2_sce_milo <- makeNhoods(num2_sce_milo, prop = 0.3, k = 30, d=8, refined = TRUE)
GSE138852_sce_milo <- makeNhoods(GSE138852_sce_milo, prop = 0.6, k = 100, d=50, refined = FALSE) # if BBKNN refined = FALSE
plotNhoodSizeHist(GSE138852_sce_milo)

##Add colData
tmp_df <- colData(GSE138852_sce_milo)
meta_file <- read.csv("./GSE138852_test/meta.disease.csv",row.names = 1)
head(meta_file)
GSE138852_meta <- meta_file[row.names(tmp_df),]
dim(GSE138852_meta)
dim(tmp_df)
head(GSE138852_meta)
GSE138852_meta_tmp <- GSE138852_meta[,c("GSE","GSM","major_cluster","disease2")]
colData(GSE138852_sce_milo) <- DataFrame(GSE138852_meta_tmp)

GSE138852_sce_milo <- countCells(GSE138852_sce_milo, meta.data = data.frame(colData(GSE138852_sce_milo)), samples="GSM")
head(nhoodCounts(GSE138852_sce_milo))


##8 differential abundance testing
GSE138852_design <- data.frame(colData(GSE138852_sce_milo))[,c("GSM","disease2")]
GSE138852_design <- distinct(GSE138852_design)
colnames(GSE138852_design) <- c("Sample","Condition")
#GSE138852_design$Condition <- c("disease","disease","control","control")
rownames(GSE138852_design) <- GSE138852_design$Sample
GSE138852_design <- GSE138852_design[colnames(nhoodCounts(GSE138852_sce_milo)), , drop=FALSE]
GSE138852_design

#GSE138852_sce_milo <- calcNhoodDistance(GSE138852_sce_milo, d=50)
rownames(GSE138852_design) <- GSE138852_design$Sample
da_results <- testNhoods(GSE138852_sce_milo, design = ~ Condition, design.df = GSE138852_design)

da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

##9 visualize neighourhoods displaying DA
GSE138852_sce_milo <- buildNhoodGraph(GSE138852_sce_milo)

plotUMAP(GSE138852_sce_milo) + plotNhoodGraphDA(GSE138852_sce_milo, da_results, alpha=0.05) +
  plot_layout(guides="collect")

##plot 
####################### 
###ann data
da_results <- annotateNhoods(GSE138852_sce_milo, da_results, coldata_col = "major_cluster")
head(da_results)

ggplot ( da_results , aes ( major_cluster_fraction )) + geom_histogram ( bins = 50 )
da_results$major_cluster <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$major_cluster)
plotDAbeeswarm(da_results, group.by = "major_cluster", alpha = 0.1)
