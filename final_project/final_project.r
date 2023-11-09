### 載入packages
library("Seurat")
library("dplyr")
library("magrittr")

# 資料path設定
setwd("/Users/benson/Desktop/Transcriptome")
data.dir <- file.path("finalproject_scRNAseq")
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
  barcode.path <- paste0(dir, "barcodes.tsv.gz")
  features.path <- paste0(dir, "features.tsv.gz")
  matrix.path <- paste0(dir, "matrix.mtx.gz")
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
# sample_data <- c("sample1","sample2","sample3","sample4","sample5","sample6")
# seurat_list <- list()
# for (i in sample_data) {
#   var_name <- paste0(i, "_seurat")
#   assign(var_name, createSeuratObjEachSample(
#     project.nam = i,
#     dir = paste0(data.dir, "/data/", i),
#     min.cells = 10,
#     min.features = 0,
#     assay = "RNA",
#     max.mt.percentage = 20,
#     num.hvf = 2000,
#     mt.pattern = "^MT-",
#     normalization.method = "LogNormalize",
#     hvf.selection.method = "vst"
#   ))
# }


sample1_seurat <- createSeuratObjEachSample (project.nam = "sample1", 
                                             dir = file.path(data.dir, "data","sample1"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^MT-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")
sample2_seurat <- createSeuratObjEachSample (project.nam = "sample2", 
                                             dir = file.path(data.dir, "data","sample2"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^MT-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")
sample3_seurat <- createSeuratObjEachSample (project.nam = "sample3", 
                                             dir = file.path(data.dir, "data","sample3"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^MT-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")
sample4_seurat <- createSeuratObjEachSample (project.nam = "sample4", 
                                             dir = file.path(data.dir, "data","sample4"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^MT-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")
sample5_seurat <- createSeuratObjEachSample (project.nam = "sample5", 
                                             dir = file.path(data.dir, "data","sample5"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^MT-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")
sample6_seurat <- createSeuratObjEachSample (project.nam = "sample6", 
                                             dir = file.path(data.dir, "data","sample6"),
                                             min.cells = 10,
                                             min.features = 0,
                                             assay = "RNA",
                                             max.mt.percentage = 20,
                                             num.hvf = 2000,
                                             mt.pattern = "^MT-",
                                             normalization.method = "LogNormalize",
                                             hvf.selection.method = "vst")

# Data integration
# Create an ‘integrated’ data assay for downstream analysis
# Identify cell types that are present in both datasets
# Obtain cell type markers that are conserved in both control and stimulated cells
# Compare the datasets to find cell-type specific responses to stimulation
data.list1 <- list('sample1' = sample1_seurat, 'sample4' = sample4_seurat)
# select features that are repeatedly variable across datasets for integration
features <- Seurat::SelectIntegrationFeatures(object.list = data.list1, verbose = FALSE)
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
DefaultAssay(combined.obj) <- "integrated"
dim(combined.obj)

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
### 建議儲存下來
saveRDS(combined.obj, "my_integrated_seurat_object.rds")
combined.obj@meta.data %>% head()
combined.obj@meta.data$orig.ident %>% table()

Seurat::DimPlot(combined.obj, reduction="tsne")
Seurat::DimPlot(combined.obj, reduction="tsne", group.by = "orig.ident")
FeaturePlot(combined.obj, reduction="tsne", features = "Mito_percent")









# # 計算每個 cluster 中細胞數量的比例
# cell_proportions <- combined.obj@meta.data[,c(1,6)]
# a <- table(cell_proportions) %>% t() %>% as.data.frame()
# control <- a$Freq[1:21] %>% sum()
# treat <- a$Freq[22:42] %>% sum()
# con <- a[1:21,]
# tre <- a[22:42,]
# con$ratio <- con$Freq/control
# tre$ratio <- tre$Freq/treat
# final <- rbind(con,tre)
# # 創建 bar plot
# library(ggplot2)
# p <- ggplot(data = final, aes(x = ratio, y = orig.ident, fill = seurat_clusters)) +
#   geom_bar(stat = "identity")
# p
# #
# markers <- FindAllMarkers(combined.obj, only.pos = TRUE)
# head(markers)
# sg <- list(cell1 = c("Rps6ka2", "Slpi", "Slfn4", "Ifitm1", "Cxcl2"),
#            cell2 = c("Mrc1", "Ly86", "Ms4a4c", "Pld4", "Dse"))
# names(combined.obj@meta.data)
# seurat.obj <- Seurat::AddModuleScore(combined.obj, features = markers$gene[1:10],
#                                      name = "sg")
# names(seurat.obj@meta.data)
# FeaturePlot(seurat.obj, reduction = "umap", features = "sg1")
# FeaturePlot(seurat.obj, reduction = "umap", features = "sg2")
