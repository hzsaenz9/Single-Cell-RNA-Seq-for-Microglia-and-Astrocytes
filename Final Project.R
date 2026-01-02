install.packages('dplyr') 
install.packages('Seurat') 
install.packages('gridExtra') 
install.packages("fs") 
install.packages('devtools') 
devtools::install_github('immunogenomics/presto') 
library(dplyr) 
library(Seurat) 
library(patchwork) 
library(gridExtra) 
library(devtools) 
library(ggplot2) 
options(java.parameters = "-Xmx8g") 
options(java.parameters = "-Xmx48g") 
# Load the PBMC dataset 
ApoE3LPS.data <- Read10X(data.dir = "C:/Users/hzepeda6/Desktop/Single Cell/APOE modulates microglial 
immunometabolism in response to age, amyloid pathology, and inflammatory 
challenge/APOExLPS/filtered_gene_bc_matrices/ApoE3LPS/hg19") 
ApoE3NaCl.data <- Read10X(data.dir = "C:/Users/hzepeda6/Desktop/Single Cell/APOE modulates microglial 
immunometabolism in response to age, amyloid pathology, and inflammatory 
challenge/APOExLPS/filtered_gene_bc_matrices/ApoE3NaCl/hg19") 
ApoE4LPS.data <- Read10X(data.dir = "C:/Users/hzepeda6/Desktop/Single Cell/APOE modulates microglial 
immunometabolism in response to age, amyloid pathology, and inflammatory 
challenge/APOExLPS/filtered_gene_bc_matrices/ApoE4LPS/hg19") 
ApoE4NaCl.data <- Read10X(data.dir = "C:/Users/hzepeda6/Desktop/Single Cell/APOE modulates microglial 
immunometabolism in response to age, amyloid pathology, and inflammatory 
challenge/APOExLPS/filtered_gene_bc_matrices/ApoE4NaCl/hg19") 
ApoE3LPS <- CreateSeuratObject(counts = ApoE3LPS.data, project = "ApoE3LPS", min.cells = 3, min.features = 
                                 200) 
ApoE3NaCl <- CreateSeuratObject(counts = ApoE3NaCl.data, project = "ApoE3NaCl", min.cells = 3, min.features 
                                = 200) 
ApoE4LPS <- CreateSeuratObject(counts = ApoE4LPS.data, project = "ApoE4LPS", min.cells = 3, min.features = 
                                 200) 
ApoE4NaCl <- CreateSeuratObject(counts = ApoE4NaCl.data, project = "ApoE4NaCl", min.cells = 3, min.features 
                                = 200) 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats 
ApoE3LPS[["percent.mt"]] <- PercentageFeatureSet(ApoE3LPS, pattern = "^MT-") 
ApoE3NaCl[["percent.mt"]] <- PercentageFeatureSet(ApoE3NaCl, pattern = "^MT-") 
ApoE4LPS[["percent.mt"]] <- PercentageFeatureSet(ApoE4LPS, pattern = "^MT-") 
ApoE4NaCl[["percent.mt"]] <- PercentageFeatureSet(ApoE4NaCl, pattern = "^MT-") 
# Visualize QC metrics as a violin plot 
#VlnPlotApoE3LPS <- VlnPlot(ApoE3LPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 
3) 
#VlnPlotApoE3NaCl <- VlnPlot(ApoE3NaCl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 
3) 
#VlnPlotApoE4LPS <- VlnPlot(ApoE4LPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 
3) 
#VlnPlotApoE4NaCl <- VlnPlot(ApoE4NaCl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 
3) 
#grid.arrange(VlnPlotApoE3LPS, VlnPlotApoE3NaCl, VlnPlotApoE4LPS, VlnPlotApoE4NaCl, nrow = 8, ncol = 
8) 
#VlnPlotApoE3LPS+VlnPlotApoE3NaCl+VlnPlotApoE4LPS+VlnPlotApoE4NaCl 
VlnPlot(ApoE3LPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(ApoE3NaCl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(ApoE4LPS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
VlnPlot(ApoE4NaCl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used 
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc. 
plot1 <- FeatureScatter(ApoE3LPS, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(ApoE3LPS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2 
ApoE3LPS <- NormalizeData(ApoE3LPS, normalization.method = "LogNormalize", scale.factor = 10000) 
ApoE3LPS <- NormalizeData(ApoE3LPS) 
ApoE3LPS <- FindVariableFeatures(ApoE3LPS, selection.method = "vst", nfeatures = 2000) 
ApoE3NaCl <- NormalizeData(ApoE3NaCl, normalization.method = "LogNormalize", scale.factor = 10000) 
ApoE3NaCl <- NormalizeData(ApoE3NaCl) 
ApoE3NaCl <- FindVariableFeatures(ApoE3NaCl, selection.method = "vst", nfeatures = 2000) 
ApoE4LPS <- NormalizeData(ApoE4LPS, normalization.method = "LogNormalize", scale.factor = 10000) 
ApoE4LPS <- NormalizeData(ApoE4LPS) 
ApoE4LPS <- FindVariableFeatures(ApoE4LPS, selection.method = "vst", nfeatures = 2000) 
ApoE4NaCl <- NormalizeData(ApoE4NaCl, normalization.method = "LogNormalize", scale.factor = 10000) 
ApoE4NaCl <- NormalizeData(ApoE4NaCl) 
ApoE4NaCl <- FindVariableFeatures(ApoE4NaCl, selection.method = "vst", nfeatures = 2000) 
# Identify the 10 most highly variable genes 
top10 <- head(VariableFeatures(ApoE3LPS), 10) 
# plot variable features with and without labels 
plot1 <- VariableFeaturePlot(ApoE3LPS) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
plot1 + plot2 
top10 <- head(VariableFeatures(ApoE3NaCl), 10) 
# plot variable features with and without labels 
plot1 <- VariableFeaturePlot(ApoE3NaCl) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
plot1 + plot2 
top10 <- head(VariableFeatures(ApoE4LPS), 10) 
# plot variable features with and without labels 
plot1 <- VariableFeaturePlot(ApoE4LPS) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
plot1 + plot2 
top10 <- head(VariableFeatures(ApoE4NaCl), 10) 
# plot variable features with and without labels 
plot1 <- VariableFeaturePlot(ApoE4NaCl) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
plot1 + plot2 
all.genes <- rownames(ApoE3LPS) 
ApoE3LPS <- ScaleData(ApoE3LPS, features = all.genes) 
ApoE3LPS <- RunPCA(ApoE3LPS, features = VariableFeatures(object = ApoE3LPS)) 
# Examine and visualize PCA results a few different ways 
print(ApoE3LPS[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(ApoE3LPS, dims = 1:2, reduction = "pca") 
DimPlot(ApoE3LPS, reduction = "pca") + NoLegend() 
DimPlotApoE3LPS <- DimPlot(ApoE3LPS, reduction = "pca") + NoLegend() 
DimHeatmap(ApoE3LPS, dims = 1:15, cells = 500, balanced = TRUE) 
#################### 
all.genes <- rownames(ApoE3NaCl) 
ApoE3NaCl <- ScaleData(ApoE3NaCl, features = all.genes) 
ApoE3NaCl <- RunPCA(ApoE3NaCl, features = VariableFeatures(object = ApoE3NaCl)) 
# Examine and visualize PCA results a few different ways 
print(ApoE3NaCl[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(ApoE3NaCl, dims = 1:2, reduction = "pca") 
DimPlot(ApoE3NaCl, reduction = "pca") + NoLegend() 
DimPlotApoE3NaCl <- DimPlot(ApoE3NaCl, reduction = "pca") + NoLegend() 
DimHeatmap(ApoE3NaCl, dims = 1:15, cells = 500, balanced = TRUE) 
################### 
all.genes <- rownames(ApoE4LPS) 
ApoE4LPS <- ScaleData(ApoE4LPS, features = all.genes) 
ApoE4LPS <- RunPCA(ApoE4LPS, features = VariableFeatures(object = ApoE4LPS)) 
# Examine and visualize PCA results a few different ways 
print(ApoE4LPS[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(ApoE4LPS, dims = 1:2, reduction = "pca") 
DimPlot(ApoE4LPS, reduction = "pca") + NoLegend() 
DimPlotApoE4LPS <- DimPlot(ApoE4LPS, reduction = "pca") + NoLegend() 
DimHeatmap(ApoE4LPS, dims = 1:15, cells = 500, balanced = TRUE) 
################## 
all.genes <- rownames(ApoE4NaCl) 
ApoE4NaCl <- ScaleData(ApoE4NaCl, features = all.genes) 
ApoE4NaCl <- RunPCA(ApoE4NaCl, features = VariableFeatures(object = ApoE4NaCl)) 
# Examine and visualize PCA results a few different ways 
print(ApoE4NaCl[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(ApoE4NaCl, dims = 1:2, reduction = "pca") 
DimPlot(ApoE4NaCl, reduction = "pca") + NoLegend() 
DimPlotApoE4NaCl <- DimPlot(ApoE4NaCl, reduction = "pca") + NoLegend() 
DimPlotApoE3LPS+DimPlotApoE3NaCl+DimPlotApoE4LPS+DimPlotApoE4NaCl 
DimHeatmap(ApoE4NaCl, dims = 1:15, cells = 500, balanced = TRUE) 
########################################################### 
ApoE3LPS <- FindNeighbors(ApoE3LPS, dims = 1:10) 
ApoE3LPS <- FindClusters(ApoE3LPS, resolution = 0.5) 
# Look at cluster IDs of the first 5 cells 
head(Idents(ApoE3LPS), 5) 
ApoE3LPS <- RunUMAP(ApoE3LPS, dims = 1:10) 
DimPlotApoE3LPS <-DimPlot(ApoE3LPS, reduction = "umap") 
DimPlotApoE3LPS 
################## 
ApoE3NaCl <- FindNeighbors(ApoE3NaCl, dims = 1:10) 
ApoE3NaCl <- FindClusters(ApoE3NaCl, resolution = 0.5) 
# Look at cluster IDs of the first 5 cells 
head(Idents(ApoE3NaCl), 5) 
ApoE3NaCl <- RunUMAP(ApoE3NaCl, dims = 1:10) 
DimPlotApoE3NaCl<-DimPlot(ApoE3NaCl, reduction = "umap") 
DimPlotApoE3NaCl 
################ 
ApoE4LPS <- FindNeighbors(ApoE4LPS, dims = 1:10) 
ApoE4LPS <- FindClusters(ApoE4LPS, resolution = 0.5) 
# Look at cluster IDs of the first 5 cells 
head(Idents(ApoE4LPS), 5) 
ApoE4LPS <- RunUMAP(ApoE4LPS, dims = 1:10) 
DimPlotApoE4LPS<-DimPlot(ApoE4LPS, reduction = "umap") 
DimPlotApoE4LPS 
################ 
ApoE4NaCl <- FindNeighbors(ApoE4NaCl, dims = 1:10) 
ApoE4NaCl <- FindClusters(ApoE4NaCl, resolution = 0.5) 
# Look at cluster IDs of the first 5 cells 
head(Idents(ApoE4NaCl), 5) 
################################### Different Groups UMAPS  
ApoE4NaCl <- RunUMAP(ApoE4NaCl, dims = 1:10) 
DimPlotApoE4NaCl <- DimPlot(ApoE4NaCl, reduction = "umap") 
DimPlotApoE4NaCl 
#DimPlotApoE3LPS+DimPlotApoE3NaCl+DimPlotApoE4LPS+DimPlotApoE4NaCl 
(DimPlotApoE3LPS+DimPlotApoE3NaCl+DimPlotApoE4LPS+DimPlotApoE4NaCl) & NoLegend() 
#grid.arrange(DimPlotApoE3LPS, DimPlotApoE3NaCl, DimPlotApoE4LPS, DimPlotApoE4NaCl, nrow = 2, ncol 
= 2,top = "UMAPs", grobs=gs) 
DimPlotApoE3LPS + plot_annotation(title = 'ApoE3+LPS') 
DimPlotApoE3NaCl + plot_annotation(title = 'ApoE3+Saline') 
DimPlotApoE4LPS + plot_annotation(title = 'ApoE4+LPS') 
DimPlotApoE4NaCl + plot_annotation(title = 'ApoE4+Saline') 
#################################### Violin Plots to visualize microglia markers in the different groups 
# Define the cluster identifiers 
cluster_ids <- c(0, 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18) 
# Initialize an empty list to store cluster markers 
cluster_markers_list <- list() 
# Iterate over each cluster identifier 
for (cluster_id in cluster_ids) { 
  # Find markers for the current cluster 
  cluster_markers <- FindMarkers(ApoE3LPS, ident.1 = cluster_id, ident.2 = NULL) 
  # Store the cluster markers in the list 
  cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] <- cluster_markers 
} 
# Print the head of cluster markers for each cluster 
for (cluster_id in cluster_ids) { 
  # Extract the cluster markers from the list 
  cluster_markers <- cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] 
  # Print the head of cluster markers 
  cat("Cluster", cluster_id, "markers:\n") 
  print(head(cluster_markers, n = 5)) 
  cat("\n") 
} 
# Define the cluster identifiers 
cluster_ids <- c(0, 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21) 
# Initialize an empty list to store cluster markers 
cluster_markers_list <- list() 
# Iterate over each cluster identifier 
for (cluster_id in cluster_ids) { 
  # Find markers for the current cluster 
  cluster_markers <- FindMarkers(ApoE3NaCl, ident.1 = cluster_id, ident.2 = NULL) 
  # Store the cluster markers in the list 
  cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] <- cluster_markers 
} 
# Print the head of cluster markers for each cluster 
for (cluster_id in cluster_ids) { 
  # Extract the cluster markers from the list 
  cluster_markers <- cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] 
  # Print the head of cluster markers 
  cat("Cluster", cluster_id, "markers:\n") 
  print(head(cluster_markers, n = 5)) 
  cat("\n") 
} 
# Define the cluster identifiers 
cluster_ids <- c(0, 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15) 
# Initialize an empty list to store cluster markers 
cluster_markers_list <- list() 
# Iterate over each cluster identifier 
for (cluster_id in cluster_ids) { 
  # Find markers for the current cluster 
  cluster_markers <- FindMarkers(ApoE4LPS, ident.1 = cluster_id, ident.2 = NULL) 
  # Store the cluster markers in the list 
  cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] <- cluster_markers 
} 
# Print the head of cluster markers for each cluster 
for (cluster_id in cluster_ids) { 
  # Extract the cluster markers from the list 
  cluster_markers <- cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] 
  # Print the head of cluster markers 
  cat("Cluster", cluster_id, "markers:\n") 
  print(head(cluster_markers, n = 5)) 
  cat("\n") 
} 
# Define the cluster identifiers 
cluster_ids <- c(0, 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19) 
# Initialize an empty list to store cluster markers 
cluster_markers_list <- list() 
# Iterate over each cluster identifier 
for (cluster_id in cluster_ids) { 
  # Find markers for the current cluster 
  cluster_markers <- FindMarkers(ApoE4NaCl, ident.1 = cluster_id, ident.2 = NULL) 
  # Store the cluster markers in the list 
  cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] <- cluster_markers 
} 
# Print the head of cluster markers for each cluster 
for (cluster_id in cluster_ids) { 
  # Extract the cluster markers from the list 
  cluster_markers <- cluster_markers_list[[paste0("cluster", cluster_id, ".markers")]] 
  # Print the head of cluster markers 
  cat("Cluster", cluster_id, "markers:\n") 
  print(head(cluster_markers, n = 5)) 
  cat("\n") 
} 
############################# ApoE3LPS 
new.cluster.idsApoe3LPS <- c("Reactive Microglia", "Oligodendrocytes", "Macrophages",  
                             "Excitatory Neurons", "Astrocytes", "Endothelial Cells", 
                             "Mature Oligodendrocytes", "Epithelial Cells", 
                             "Pyramidal Neurons", "Microglia","Pericytes","Radial Glia", 
                             "Vascular Smooth Muscle Cells","Fibroblasts", "Reactive Astrocytes", 
                             "Reactive Astrocytes","Migrating Neuroblasts", "Reactive Microglia",  
                             "Endothelial Cells") 
# Assign cluster names to Seurat object 
names(new.cluster.idsApoe3LPS) <- levels(ApoE3LPS) 
# Rename cluster identities 
ApoE3LPS <- RenameIdents(ApoE3LPS, new.cluster.idsApoe3LPS) 
# Plot UMAP with cluster names 
DimPlot(ApoE3LPS, reduction = "umap", label = TRUE, pt.size = 0.5) +  plot_annotation(title = 'ApoE3+LPS') 
# Visualize gene expression using FeaturePlot 
FeaturePlot(ApoE3LPS, features = c("Serping1", "Gfap")) +  plot_annotation(title = 'ApoE3+LPS') 
FeaturePlot(ApoE3LPS, features = c("Ctsb","Ctsd","Ctsl")) +  plot_annotation(title = 'ApoE3+LPS') 
FeaturePlot(ApoE3LPS, features = c("Plaur","Mmp2","Mmp13")) +  plot_annotation(title = 'ApoE3+LPS') 
FeaturePlot(ApoE3LPS, features = c("Axin2","Nes","Ctnnb1")) +  plot_annotation(title = 'ApoE3+LPS') 
VlnPlot(ApoE3LPS, features = c("Serping1", "Gfap")) +  plot_annotation(title = 'ApoE3+LPS') 
VlnPlot(ApoE3LPS, features = c("Ctsb", "Ctsl")) +  plot_annotation(title = 'ApoE3+LPS') 
VlnPlot(ApoE3LPS, features = c("Trem2","Tyrobp")) +  plot_annotation(title = 'ApoE3+LPS') 
VlnPlot(ApoE3LPS, features = c("Aebp1","Wwtr1")) +  plot_annotation(title = 'ApoE3+LPS') 
VlnPlot(ApoE3LPS, features = c("Phyhd1","Dst","Rasl12")) +  plot_annotation(title = 'ApoE3+LPS') 
VlnPlot(ApoE3LPS, features = c("Olig1")) +  plot_annotation(title = 'ApoE3+LPS') 
############################################Apoe3+Saline 
new.cluster.idsApoE3NaCl <- c("Astrocytes", "Excitatory Neurons", "Microglia",  
                              "Mature Oligodendrocytes", "Oligodendrocytes", "Endothelial Cells", 
                              "Oligodendrocytes", "Microglia", "Vascular Smooth Muscle Cells", 
                              "Vascular Smooth Muscle Cells","Neuroblasts","Endothelial Cells", 
                              "Astrocytes","Microglia","Reactive Microglia","Mature Oligodendrocytes", 
                              "Endothelial Cells","Mature Oligodendrocytes","Astrocytes","Fibroblasts", 
                              "Neuroblasts","Endothelial Cells") 
# Assign cluster names to Seurat object 
names(new.cluster.idsApoE3NaCl) <- levels(ApoE3NaCl) 
# Rename cluster identities 
ApoE3NaCl <- RenameIdents(ApoE3NaCl, new.cluster.idsApoE3NaCl) 
# Plot UMAP with cluster names 
DimPlot(ApoE3NaCl, reduction = "umap", label = TRUE, pt.size = 0.5) +  plot_annotation(title = 'ApoE3+Saline') 
# Visualize gene expression using FeaturePlot 
FeaturePlot(ApoE3NaCl, features = c("Serping1", "Gfap", "Ggta")) +  plot_annotation(title = 'ApoE3+Saline') 
FeaturePlot(ApoE3NaCl, features = c("Ctsb","Ctsd","Ctsl")) +  plot_annotation(title = 'ApoE3+Saline') 
FeaturePlot(ApoE3NaCl, features = c("Plaur","Mmp2","Mmp13")) +  plot_annotation(title = 'ApoE3+Saline') 
FeaturePlot(ApoE3NaCl, features = c("Axin2","Nes","Ctnnb1")) +  plot_annotation(title = 'ApoE3+Saline') 
# Get the gene expression matrix 
expression_matrix <- GetAssayData(ApoE3NaCl, assay = "RNA") 
# Extract gene expression values for specific genes 
gene_expression <- expression_matrix[c("Serping1", "Gfap"), ] 
VlnPlot(ApoE3NaCl, features = c("Serping1", "Gfap")) +  plot_annotation(title = 'ApoE3+Saline') 
VlnPlot(ApoE3NaCl, features = c("Ctbs", "Ctsl")) +  plot_annotation(title = 'ApoE3+Saline') 
VlnPlot(ApoE3NaCl, features = c("Trem2","Tyrobp")) +  plot_annotation(title = 'ApoE3+Saline') 
VlnPlot(ApoE3NaCl, features = c("Aebp1","Wwtr1")) +  plot_annotation(title = 'ApoE3+Saline') 
VlnPlot(ApoE3NaCl, features = c("Phyhd1","Dst","Rasl12")) +  plot_annotation(title = 'ApoE3+Saline') 
VlnPlot(ApoE3NaCl, features = c("Olig1")) +  plot_annotation(title = 'ApoE3+Saline') 
#####################################Apoe4+LPS 
new.cluster.idsApoE4LPS <- c("Glutamergic Neurons","Microglia","Oligodendrocytes", 
                             "Reactive Astrocytes","Oligodendrocytes","Endothelial Cells", 
                             "Neurons","Neurons","Endothelial Cells", 
                             "Astrcoytes","Neurons","Oligodendrocytes", 
                             "Oligodendrocytes","Oligodendrocytes","Microglia", 
                             "Reactive Astrocytes","Neural Progenitor Cells","Monocytes", 
                             "Endothelial Cells","Endothelial Cells") 
# Assign cluster names to Seurat object 
names(new.cluster.idsApoE4LPS) <- levels(ApoE4LPS) 
# Rename cluster identities 
ApoE4LPS <- RenameIdents(ApoE4LPS, new.cluster.idsApoE4LPS) 
# Plot UMAP with cluster names 
DimPlot(ApoE4LPS, reduction = "umap", label = TRUE, pt.size = 0.5) +  plot_annotation(title = 'ApoE4+LPS') 
# Visualize gene expression using FeaturePlot 
FeaturePlot(ApoE4LPS, features = c("Serping1", "Gfap", "Ggta")) +  plot_annotation(title = 'ApoE4+LPS') 
FeaturePlot(ApoE4LPS, features = c("Ctsb","Ctsd","Ctsl")) +  plot_annotation(title = 'ApoE4+LPS') 
FeaturePlot(ApoE4LPS, features = c("Plaur","Mmp2","Mmp13")) +  plot_annotation(title = 'ApoE4+LPS') 
FeaturePlot(ApoE4LPS, features = c("Axin2","Nes","Ctnnb1")) +  plot_annotation(title = 'ApoE4+LPS') 
VlnPlot(ApoE4LPS, features = c("Serping1", "Gfap")) +  plot_annotation(title = 'ApoE4+LPS') 
VlnPlot(ApoE4LPS, features = c("Ctsb", "Ctsl")) +  plot_annotation(title = 'ApoE4+LPS') 
VlnPlot(ApoE4LPS, features = c("Trem2","Tyrobp")) +  plot_annotation(title = 'ApoE4+LPS') 
VlnPlot(ApoE4LPS, features = c("Aebp1","Wwtr1")) +  plot_annotation(title = 'ApoE4+LPS') 
VlnPlot(ApoE4LPS, features = c("Phyhd1","Dst","Rasl12")) +  plot_annotation(title = 'ApoE4+LPS') 
VlnPlot(ApoE4LPS, features = c("Olig1")) +  plot_annotation(title = 'ApoE4+LPS') 
##################################Apoe4+Saline 
new.cluster.idsApoE4NaCl <- c("Glutamergic Neurons", "Reactive Microglia", "Microglia",  
                              "Microglia", "Oligodendrocytes", "Endothelial Cells",  
                              "Astrocytes", "Neural Progenitor Cells", "Endothelial Cells", 
                              "Fibroblasts", "Endothelial Cells", "Reactive Astrocytes", 
                              "Oligodendrocytes", "Oligodendrocytes", "Macrophages", 
                              "Fibroblasts", "Neural Progenitor Cells", "Endothelial Cells", 
                              "Fibroblasts", "Endothelial Cells") 
# Assign cluster names to Seurat object 
names(new.cluster.idsApoE4NaCl) <- levels(ApoE4NaCl) 
# Rename cluster identities 
ApoE4NaCl <- RenameIdents(ApoE4NaCl, new.cluster.idsApoE4NaCl) 
# Plot UMAP with cluster names 
DimPlot(ApoE4NaCl, reduction = "umap", label = TRUE, pt.size = 0.5) +  plot_annotation(title = 'ApoE4+Saline') 
FeaturePlot(ApoE4NaCl, features = c("Serping1", "Gfap", "Ggta"))  +  plot_annotation(title = 'ApoE4+Saline') 
FeaturePlot(ApoE4NaCl, features = c("Ctsb","Ctsd","Ctsl")) +  plot_annotation(title = 'ApoE4+NaCl') 
FeaturePlot(ApoE4NaCl, features = c("Plaur","Mmp2","Mmp13")) +  plot_annotation(title = 'ApoE4+Saline') 
FeaturePlot(ApoE4NaCl, features = c("Axin2","Nes","Ctnnb1")) +  plot_annotation(title = 'ApoE4+Saline') 
VlnPlot(ApoE4NaCl, features = c("Serping1", "Gfap")) +  plot_annotation(title = 'ApoE4+Saline') 
VlnPlot(ApoE4NaCl, features = c("Ctsb", "Ctsl")) +  plot_annotation(title = 'ApoE4+Saline') 
VlnPlot(ApoE4NaCl, features = c("Trem2","Tyrobp")) +  plot_annotation(title = 'ApoE4+Saline') 
VlnPlot(ApoE4NaCl, features = c("Aebp1","Wwtr1")) +  plot_annotation(title = 'ApoE4+Saline') 
VlnPlot(ApoE4NaCl, features = c("Phyhd1","Dst","Rasl12")) +  plot_annotation(title = 'ApoE4+Saline') 
VlnPlot(ApoE4NaCl, features = c("Olig1")) +  plot_annotation(title = 'ApoE4+Saline') 
############################################### 
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE);  
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE) 
RidgePlot(ApoE4LPS, features = c("Gfap"), ncol = 1)+  plot_annotation(title = 'ApoE4+LPS') 
RidgePlot(ApoE3LPS, features = c("Gfap"), ncol = 1)+  plot_annotation(title = 'ApoE3+LPS') 
RidgePlot(ApoE3NaCl, features = c("Gfap"), ncol = 1)+  plot_annotation(title = 'ApoE3+Saline') 
RidgePlot(ApoE4NaCl, features = c("Gfap"), ncol = 1)+  plot_annotation(title = 'ApoE4+Saline') 
RidgePlot(ApoE4LPS, features = c("Wwtr1"), ncol = 1)+  plot_annotation(title = 'ApoE4+LPS') 
RidgePlot(ApoE3LPS, features = c("Wwtr1"), ncol = 1)+  plot_annotation(title = 'ApoE3+LPS') 
RidgePlot(ApoE3NaCl, features = c("Wwtr1"), ncol = 1)+  plot_annotation(title = 'ApoE3+Saline') 
RidgePlot(ApoE4NaCl, features = c("Wwtr1"), ncol = 1)+  plot_annotation(title = 'ApoE4+Saline') 
RidgePlot(ApoE4LPS, features = c("Trem2"), ncol = 1)+  plot_annotation(title = 'ApoE4+LPS') 
RidgePlot(ApoE3LPS, features = c("Trem2"), ncol = 1)+  plot_annotation(title = 'ApoE3+LPS') 
RidgePlot(ApoE3NaCl, features = c("Trem2"), ncol = 1)+  plot_annotation(title = 'ApoE3+Saline') 
RidgePlot(ApoE4NaCl, features = c("Trem2"), ncol = 1)+  plot_annotation(title = 'ApoE4+Saline') 
RidgePlot(ApoE4LPS, features = c("Ctsb"), ncol = 2) 
#################################### Combined Groups 
pbmc.combinedApoE3 <- merge(ApoE3NaCl, y = ApoE3LPS, add.cell.ids = c("15251", "15927"), project = 
                              "PBMC12K") 
pbmc.combinedApoE3 
pbmc.combinedApoE4 <- merge(ApoE4NaCl, y = ApoE4LPS, add.cell.ids = c("13911", "29259"), project = 
                              "PBMC12K") 
pbmc.combinedApoE4 
pbmc.combinedApoE <- merge(pbmc.combinedApoE3, y = pbmc.combinedApoE4, add.cell.ids = c("31178", 
                                                                                        "43170"), project = "PBMC12K") 
pbmc.combinedApoE 
################################################### 
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats 
pbmc.combinedApoE[["percent.mt"]] <- PercentageFeatureSet(pbmc.combinedApoE, pattern = "^MT-") 
# Visualize QC metrics as a violin plot 
VlnPlot(pbmc.combinedApoE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used 
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc. 
plot1 <- FeatureScatter(pbmc.combinedApoE, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(pbmc.combinedApoE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2 
pbmc.combinedApoE <- NormalizeData(pbmc.combinedApoE, normalization.method = "LogNormalize", 
                                   scale.factor = 10000) 
pbmc.combinedApoE <- NormalizeData(pbmc.combinedApoE) 
pbmc.combinedApoE <- FindVariableFeatures(pbmc.combinedApoE, selection.method = "vst", nfeatures = 2000) 
# Identify the 10 most highly variable genesz 
top10 <- head(VariableFeatures(pbmc.combinedApoE), 10) 
# plot variable features with and without labels 
plot1 <- VariableFeaturePlot(pbmc.combinedApoE) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) 
plot1 + plot2 
all.genes <- rownames(pbmc.combinedApoE) 
pbmc.combinedApoE <- ScaleData(pbmc.combinedApoE, features = all.genes) 
pbmc.combinedApoE <- RunPCA(pbmc.combinedApoE, features = VariableFeatures(object = 
                                                                             pbmc.combinedApoE)) 
# Examine and visualize PCA results a few different ways 
print(pbmc.combinedApoE[["pca"]], dims = 1:5, nfeatures = 5) 
VizDimLoadings(pbmc.combinedApoE, dims = 1:2, reduction = "pca") 
DimPlot(pbmc.combinedApoE, reduction = "pca") + NoLegend() 
DimHeatmap(pbmc.combinedApoE, dims = 1:15, cells = 500, balanced = TRUE) 
pbmc.combinedApoE <- FindNeighbors(pbmc.combinedApoE, dims = 1:10) 
pbmc.combinedApoE <- FindClusters(pbmc.combinedApoE, resolution = 0.5) 
# Look at cluster IDs of the first 5 cells 
head(Idents(pbmc.combinedApoE), 5) 
pbmc.combinedApoE <- RunUMAP(pbmc.combinedApoE, dims = 1:10) 
DimPlot(pbmc.combinedApoE, reduction = "umap") 
pbmc.combinedApoE <- JoinLayers(pbmc.combinedApoE) 
RidgePlot(pbmc.combinedApoE, features = c("Gfap"), ncol = 1)+  plot_annotation(title = 'Combined ApoE') 
RidgePlot(pbmc.combinedApoE, features = c("Wwtr1"), ncol = 1)+  plot_annotation(title = 'Combined ApoE') 
RidgePlot(pbmc.combinedApoE, features = c("Aebp1"), ncol = 1)+  plot_annotation(title = 'Combined ApoE') 
RidgePlot(pbmc.combinedApoE, features = c("Trem2"), ncol = 1)+  plot_annotation(title = 'Combined ApoE') 
file.copy(from=plots.png.paths, to="C:/Users/hzepeda6/Desktop/Single Cell/APOE modulates microglial 
immunometabolism in response to age, amyloid pathology, and inflammatory challenge/New New Plots") 