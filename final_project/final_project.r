### 載入packages
library("Seurat")
library("dplyr")
library("magrittr")

# 資料path設定
setwd("/Users/benson/Desktop/Transcriptome")
data.dir <- file.path("final")
# 將單一樣本資料前處理，寫成一個function
# 讀檔 (matrix.mtx, barcode.tsv, and feature.tsv)
# 去除duplicated genes
# 資料載入seurat
# 計算mito gene的比例
# Cell QC
# Normalization
# Identify highly variable genes
createSeuratObjEachSample <- function(project.name, 
                                      dir,
                                      min.cells = 10,
                                      min.features = 0,
                                      assay = "RNA",
                                      max.mt.percentage = 20,
                                      num.hvf = 2000,
                                      mt.pattern = "^mt-",
                                      normalization.method = "LogNormalize",
                                      hvf.selection.method = "vst"){
  #### 1. Read files
  cat(c("Read files\n"))
  barcode.path <- file.path(dir, "barcodes.tsv")
  features.path <- file.path(dir, "genes.tsv")
  matrix.path <- file.path(dir, "matrix.mtx")
  mat <- Matrix::readMM(file = matrix.path)
  feature.names = read.delim(features.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  barcode.names = read.delim(barcode.path,
                             header = FALSE,
                             stringsAsFactors = FALSE)
  colnames(mat) = paste(project.name, barcode.names$V1, sep = "_")
  rownames(mat) = feature.names$V2
  mat <- as.matrix(mat)
  cat(c("-> matrix:", dim(mat), "\n"))
  
  #### 2. 將基因名子一樣的資料加總起來
  cat(c("Remove duplicated genes\n"))
  IDfreqs <- table(rownames(mat))
  IDfreqs <- IDfreqs[IDfreqs > 1]
  cat(c("-> Num. dup. genes:", length(IDfreqs), "\n"))
  if (length(IDfreqs) > 0){
    for (i in names(IDfreqs)) {
      index <- which(rownames(mat) == i)
      duplicates <- as.matrix(mat[index, ])
      sums <- apply(duplicates, 2, sum, na.rm=TRUE)
      
      ### 置換成加總數值
      mat[index[1],] <- sums
      
      ### 將重複的資料移除
      for (j in length(index):2) {
        mat <- mat[-index[j], ]
      }
    }
  }
  cat(c("-> matrix:", dim(mat), "\n"))
  
  #### 3. 資料載入seurat
  cat(c("Create Seurat Object\n"))
  seurat.obj <- Seurat::CreateSeuratObject(counts = mat, 
                                           project = project.name,
                                           min.cells = min.cells, ## Include features detected in at least this many cells.
                                           min.features = 0, ## Include cells where at least this many features are detected.
                                           assay = assay, ## Name of the assay corresponding to the initial input data.
                                           meta.data = NULL)
  cat(c("-> dimension:", dim(seurat.obj), "\n"))
  
  #### 4. 計算mito gene的比例
  cat(c("Calculate MT percentage\n"))
  seurat.obj[["Mito_percent"]] <- Seurat::PercentageFeatureSet(seurat.obj, pattern = mt.pattern)
  
  ##### 5. Cell QC
  cat(c("Cell QC\n"))
  min_nFeature <- min.features
  max_Mito_percent <- max.mt.percentage
  seurat.obj <- seurat.obj %>%
    subset(., subset = nFeature_RNA > min_nFeature) %>%
    subset(., subset = Mito_percent < max_Mito_percent)
  cat(c("-> dimension:", dim(seurat.obj), "\n"))
  
  ##### 6. Normalization
  cat(c("Normalization\n"))
  seurat.obj <- Seurat::NormalizeData(seurat.obj,
                                      normalization.method = normalization.method)
  
  #### 7. Identify highly variable genes
  cat(c("High variable genes\n"))
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, 
                                             selection.method = hvf.selection.method, 
                                             nfeatures = num.hvf #Number of features to select as top variable features
  )
  
  return(seurat.obj)
}

# 1. 將每個樣本建構出Seurat object
k_seurat <- createSeuratObjEachSample (project.nam = "k", 
                                             dir = file.path(data.dir, "K"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^mt-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")
kl_seurat <- createSeuratObjEachSample (project.nam = "kl", 
                                             dir = file.path(data.dir, "KL"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^mt-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")


# Data integration
# Create an ‘integrated’ data assay for downstream analysis
# Identify cell types that are present in both datasets
# Obtain cell type markers that are conserved in both control and stimulated cells
# Compare the datasets to find cell-type specific responses to stimulation
data.list <- list("k" = k_seurat, "kl" = kl_seurat)
# select features that are repeatedly variable across datasets for integration
features <- Seurat::SelectIntegrationFeatures(object.list = data.list, 
                                              nfeatures = 1000,
                                              verbose = FALSE)
class(features)
length(features)
features[1:10]
#### identify anchors 
anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, 
                                          anchor.features = features,
                                          verbose = F)


# this command creates an 'integrated' data assay
combined.obj <- Seurat::IntegrateData(anchorset = anchors)
DefaultAssay(combined.obj) <- "integrated"
dim(combined.obj)

# Run the standard workflow for visualization and clustering
combined.obj <- Seurat::ScaleData(combined.obj, verbose = FALSE)
combined.obj <- Seurat::RunPCA(combined.obj, npcs = 50, verbose = FALSE)
ElbowPlot(combined.obj, ndims = 50)
pcs <- 20
combined.obj <- combined.obj %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:pcs, verbose = FALSE) %>%
  Seurat::RunTSNE(reduction="pca", dims= 1:pcs, verbose = FALSE) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:pcs, verbose = FALSE) %>%
  Seurat::FindClusters(resolution = 0.5, verbose = FALSE)
  Seurat::FindClusters(resolution = 1, verbose = FALSE)
### 建議儲存下來
saveRDS(combined.obj, "integrated_seurat_object.rds")

###### 從此開始讀檔分析～ ------------------------------------------------------
### 讀檔 
combined.obj <- file.path(data.dir, "integrated_seurat_object.rds") %>% readRDS()
combined.obj@meta.data %>% head()
combined.obj@meta.data$orig.ident %>% table()

Seurat::DimPlot(combined.obj, reduction="tsne")
Seurat::DimPlot(combined.obj, reduction="umap")
Seurat::DimPlot(combined.obj, reduction="tsne", group.by = "orig.ident")
FeaturePlot(combined.obj, reduction="tsne", features = "Mito_percent")


# # Assigning cell type identity to clusters
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet","10","11","12","13","14","15","16","17","18","19")
# names(new.cluster.ids) <- levels(combined.obj)
# com <- RenameIdents(combined.obj, new.cluster.ids)
# DimPlot(com, reduction = "umap", label = TRUE, pt.size = 0.5)


combined.obj[[]] %>% head()
Seurat::Idents(combined.obj)[1:10]
Seurat::Idents(combined.obj) <- "integrated_snn_res.0.5"
Seurat::DimPlot(combined.obj, reduction="umap", group.by = "integrated_snn_res.0.5")
Seurat::Idents(combined.obj) <- "integrated_snn_res.1"
Seurat::DimPlot(combined.obj, reduction="umap", group.by = "integrated_snn_res.1")

# # 手動更改cluster
# combined.obj@meta.data$new_column <- combined.obj@meta.data$RNA_snn_res.0.5
# combined.obj@meta.data$new_column[combined.obj@meta.data$new_column == 5] <- 0
# table(combined.obj@meta.data$new_column)

# 13. Cell doublet prediction -------------------------------------------
# 載入scDblFinder
library("scDblFinder")
count.df <- Seurat::GetAssayData(combined.obj, layer = "counts", assay = "RNA")
dbl <- scDblFinder(count.df, 
                   clusters = Seurat::Idents(combined.obj),
                   verbose = FALSE)

combined.obj[["dbl.class"]] <- dbl$scDblFinder.class
combined.obj[["dbl.score"]] <- dbl$scDblFinder.score
head(dbl)
combined.obj[[]] %>% head()
DimPlot(combined.obj, reduction="tsne", group.by = "dbl.class")

# 14. Cell cycle phase prediction and regression ----------------------------
cc.genes  # 內建cell cycle gene list
#
mouse.cc.genes <- list(
  s.genes = stringr::str_to_title(cc.genes$s.genes),
  g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)
)
mouse.cc.genes
#### Predict cell cycle phase
combined.obj <- Seurat::CellCycleScoring(combined.obj, 
                                       s.features = mouse.cc.genes, 
                                       g2m.features = mouse.cc.genes, 
                                       set.ident = FALSE)
# 觀察是否需要Regression
Seurat::DimPlot(combined.obj, reduction="umap", group.by = "Phase")

# Regression out of cell cycle phase
regress.obj <- Seurat::ScaleData(combined.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(combined.obj)) %>%
  Seurat::RunPCA(npcs = 20, verbose = FALSE) %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
  Seurat::RunTSNE(reduction="pca", dims= 1:20, verbose = FALSE) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:pcs, verbose = FALSE) %>%
  Seurat::FindClusters(resolution = 0.5, verbose = FALSE)
Seurat::DimPlot(regress.obj, reduction="umap", group.by = "Phase")

# 13. Identify marker genes & cluster annotation
# suggested reading: https://www.nature.com/articles/s41596-021-00534-0
combined.obj <- readRDS("combined.obj.RDS")
DimPlot(combined.obj)
### Cluster 2 vs. the others
c2.markers <- Seurat::FindMarkers(combined.obj, 
                                  ident.1 = 2, min.pct = 0.25,  # 至少要多少比例細胞符合
                                  only.pos = TRUE) ### Using default assay
head(c2.markers)
VlnPlot(combined.obj, features = c("Cd177"))
FeaturePlot(combined.obj, features = c("Cd177"), reduction = "umap")
# Markers for each cluster
markers <- FindAllMarkers(combined.obj, only.pos = TRUE) %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>% 
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(combined.obj, features = top10$gene) + NoLegend()
head(markers)
write.table(top10[,6:7],file = "top10.csv",sep=",",row.names = F)

# Compare c0 vs. c3
deg_bw_c0_and_c3 <- Seurat::FindMarkers(combined.obj, 
                                        ident.1 = 0, 
                                        ident.2 = 3,
                                        min.pct = 0.25, 
                                        only.pos = FALSE) ### Using default assay
head(deg_bw_c0_and_c3)
VlnPlot(combined.obj, features = c("Fbn1"))
FeaturePlot(combined.obj, features = c("Fbn1"), reduction = "umap")
# Cell marker database
# https://panglaodb.se/
#   http://bio-bigdata.hrbmu.edu.cn/CellMarker/
#   http://cloud.capitalbiotech.com/SingleCellBase/
#   Using a set of genes (signature)
sg <- list(cell1 = c("S100a8", "S100a9", "Cd177", "Mmp9", "Itgam"),
           cell2 = c("Mrc1", "Ly86", "Ms4a4c", "Pld4", "Dse"))
names(combined.obj@meta.data)
combined.obj <- Seurat::AddModuleScore(combined.obj, features = sg,
                                     name = "sg")
names(combined.obj@meta.data)
FeaturePlot(combined.obj, reduction = "umap", features = "sg1")
FeaturePlot(combined.obj, reduction = "umap", features = "sg2")
