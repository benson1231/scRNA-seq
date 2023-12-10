### 載入packages ---------------------------------------------------------------
library("Seurat")
library("dplyr")
library("magrittr")

### 將單一樣本資料前處理，寫成一個function -------------------------------------
# (1)讀檔 (matrix.mtx, barcode.tsv, and feature.tsv)
# (2)去除duplicated genes
# (3)資料載入seurat
# (4)計算mito gene的比例
# (5)Cell QC
# (6)Normalization
# (7)Identify highly variable genes
createSeuratObjEachSample <- function(project.name, 
                                      dir,
                                      min.cells = 3,
                                      min.features = 500,
                                      max.features = 6000,
                                      assay = "RNA",
                                      max.mt.percentage = 20,
                                      num.hvf = 2000,
                                      mt.pattern = "^mt-",
                                      normalization.method = "LogNormalize",
                                      hvf.selection.method = "vst"){
  #### (1) Read files
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
  
  #### (2) 將基因名子一樣的資料加總起來
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
  
  #### (3) 資料載入seurat
  cat(c("Create Seurat Object\n"))
  seurat.obj <- Seurat::CreateSeuratObject(counts = mat, 
                                           project = project.name,
                                           min.cells = min.cells, ## Include features detected in at least this many cells.
                                           min.features = min.features, ## Include cells where at least this many features are detected.
                                           assay = assay, ## Name of the assay corresponding to the initial input data.
                                           meta.data = NULL)
  cat(c("-> dimension:", dim(seurat.obj), "\n"))
  
  #### (4) 計算mito gene的比例
  cat(c("Calculate MT percentage\n"))
  seurat.obj[["Mito_percent"]] <- Seurat::PercentageFeatureSet(seurat.obj, pattern = mt.pattern)
  
  #### (5) Cell QC
  cat(c("Cell QC\n"))
  min_nFeature <- min.features
  max_nFeature <- max.features
  max_Mito_percent <- max.mt.percentage
  seurat.obj <- seurat.obj %>%
    subset(., subset = nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature) %>%
    subset(., subset = Mito_percent < max_Mito_percent)
  cat(c("-> dimension:", dim(seurat.obj), "\n"))
  
  #### (6) Normalization
  cat(c("Normalization\n"))
  seurat.obj <- Seurat::NormalizeData(seurat.obj,
                                      normalization.method = normalization.method,
                                      scale.factor = 10000)
  
  #### (7) Identify highly variable genes
  cat(c("High variable genes\n"))
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, 
                                             selection.method = hvf.selection.method, 
                                             nfeatures = num.hvf #Number of features to select as top variable features
  )
  
  return(seurat.obj)
}

### 1. 將每個樣本建構出Seurat object -------------------------------------------
k_seurat <- createSeuratObjEachSample (project.nam = "k", 
                                       dir = "K",
                                       min.cells = 3,
                                       min.features = 500,
                                       max.features = 6000,
                                       assay = "RNA",
                                       max.mt.percentage = 20,
                                       num.hvf = 2000,
                                       mt.pattern = "^mt-",
                                       normalization.method = "LogNormalize",
                                       hvf.selection.method = "vst")
kl_seurat <- createSeuratObjEachSample (project.nam = "kl", 
                                        dir = "KL",
                                        min.cells = 3,
                                        min.features = 500,
                                        max.features = 6000,
                                        assay = "RNA",
                                        max.mt.percentage = 20,
                                        num.hvf = 2000,
                                        mt.pattern = "^mt-",
                                        normalization.method = "LogNormalize",
                                        hvf.selection.method = "vst")
# 觀察兩個Seurat object
VlnPlot(k_seurat, features = c("nFeature_RNA", "nCount_RNA", "Mito_percent"), ncol = 3)
VlnPlot(kl_seurat, features = c("nFeature_RNA", "nCount_RNA", "Mito_percent"), ncol = 3)

### 2. Data integration --------------------------------------------------------
# Create an ‘integrated’ data assay for downstream analysis
# Identify cell types that are present in both datasets
# Obtain cell type markers that are conserved in both control and stimulated cells
# Compare the datasets to find cell-type specific responses to stimulation
data.list <- list("k" = k_seurat, "kl" = kl_seurat)
# select features that are repeatedly variable across datasets for integration
features <- Seurat::SelectIntegrationFeatures(object.list = data.list, 
                                              nfeatures = 2000,
                                              verbose = FALSE)
# identify anchors 
anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, 
                                          anchor.features = features,
                                          verbose = F)
# this command creates an 'integrated' data assay
combined.obj <- Seurat::IntegrateData(anchorset = anchors)
DefaultAssay(combined.obj) <- "integrated"
dim(combined.obj)
combined.obj@meta.data$orig.ident %>% table()
combined.obj@meta.data %>% head()

### 3. Run the standard workflow for visualization and clustering --------------
combined.obj <- Seurat::ScaleData(combined.obj, verbose = FALSE)
combined.obj <- Seurat::RunPCA(combined.obj, npcs = 50, verbose = FALSE)
ElbowPlot(combined.obj, ndims = 50)
pcs <- 20   # 肉眼判斷為20(符合paper作法)
# visualizing both cells and features that define the PCA
VizDimLoadings(combined.obj, dims = 1:2, reduction = "pca")
DimPlot(combined.obj, reduction = "pca") + NoLegend()
DimHeatmap(combined.obj, dims = 1, cells = 500, balanced = TRUE)

### 4. Run non-linear dimensional reduction (UMAP/tSNE) ------------------------
combined.obj <- combined.obj %>%
  Seurat::RunUMAP(reduction = "pca", umap.method = "uwot", dims = 1:pcs, verbose = FALSE) %>%
  Seurat::RunTSNE(reduction="pca", dims= 1:pcs, verbose = FALSE) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:pcs, verbose = FALSE) %>%
  Seurat::FindClusters(resolution = 0.06, verbose = FALSE) %>%  # paper提供resolution = 0.06
  Seurat::FindSubCluster("0",resolution = 0.06, graph.name = "integrated_snn", 
                         subcluster.name = "new_cluster")  # T cell進一步subcluster
# 觀察cluster狀況
Seurat::DimPlot(combined.obj, reduction="umap")  
Seurat::DimPlot(combined.obj, reduction="tsne")
Seurat::DimPlot(combined.obj, reduction="tsne", group.by = "orig.ident")
Seurat::DimPlot(combined.obj, reduction="umap", group.by = "new_cluster")
# 觀察Feature狀況
FeaturePlot(combined.obj, reduction="umap", features = "Mito_percent")

### 5. Find markers for each cluster -------------------------------------------
markers <- FindAllMarkers(combined.obj, only.pos = TRUE) %>% 
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>% 
  dplyr::arrange(desc(avg_log2FC)) %>% 
  slice_head(n = 10) %>%
  ungroup() -> top10
# 將找到的marker進行視覺化，作heatmap
Seurat::DoHeatmap(combined.obj, group.by = "ident", features = top10$gene, 
                  draw.lines = TRUE) + NoLegend()
head(markers)
# 可以將特定gene作圖
FeaturePlot(combined.obj, reduction="umap", features = "Cd19", label = T)
VlnPlot(combined.obj, features = c("Cd19", "Cd3e"))

### 6.annotation(from paper) -----------------------------------------------------
b_cell_mark <- c("Cd19","Cd22","Cd79a","Cd79b")
cancer_mark <- c("Epcam","Krt18","Krt8","Krt19")
dendritic_mark <-  c("Cd86","Cd83","Clec10a","Ccl22","Cxcl16","Xcr1","Rab7b","Naaa","Btla")
endo_mark <-  c("Pecam1","Tek","Plvap","Thbd","Tmem100","Adgrf5","Podxl","Nostrin","Acvrl1")
macrophage_mark <- c("Adgre1","Itgam","Itgal","Cd14","Ccl9","F13a1","Cx3cr1","Cd68")
neutrophil_mark <- c("Csf3r","Cxcr2S100a8")
nk_mark <- c("Klrd1","Nkg7","Prf1","Gzma","Gzmb","Ccl5","Cma1","Klre1","Klra4")
plasma_dend_mark <- c("Ccr9","Bst2")
t_cell_mark <- c("Cd3d","Cd5","Cd3e")
# 觀察cell cluster marker是否能夠代表該cluster
FeaturePlot(combined.obj, reduction="umap", features = b_cell_mark[1:4], label = T)
FeaturePlot(combined.obj, reduction="umap", features = cancer_mark[1:4], label = T)
FeaturePlot(combined.obj, reduction="umap", features = dendritic_mark[1:4], label = T)
FeaturePlot(combined.obj, reduction="umap", features = endo_mark[1:4], label = T)
FeaturePlot(combined.obj, reduction="umap", features = macrophage_mark[1:4], label = T)
FeaturePlot(combined.obj, reduction="umap", features = neutrophil_mark, label = T)
FeaturePlot(combined.obj, reduction="umap", features = nk_mark[1:4], label = T)
FeaturePlot(combined.obj, reduction="umap", features = plasma_dend_mark, label = T)
# T cell
p1 <- FeaturePlot(combined.obj, reduction="umap", features = "Cd3e", label = T)  # T cell
p2 <- FeaturePlot(combined.obj, reduction="umap", features = "Tcf7", label = T)  # activated T cell
p3 <- FeaturePlot(combined.obj, reduction="umap", features = "Ctla4", label = T) # exhausted T cell
cowplot::plot_grid(p1, p2, p3, ncol = 3)
VlnPlot(combined.obj, features = c("Cd3e","Tcf7","Ctla4"))

# 將marker_list輸出
marker_list <- c(b_cell_mark, cancer_mark, dendritic_mark, endo_mark, macrophage_mark, 
                 neutrophil_mark, nk_mark, plasma_dend_mark, t_cell_mark)
for(i in marker_list){
  FeaturePlot(combined.obj, reduction="umap", features = i, label = T) %>% print()
  ggplot2::ggsave(paste0("plot_", gsub(" ", "_", i), ".png"))
}
# cluster annotation
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "0_0"] <- "T cell activated"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "0_1"] <- "T cell exhausted"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "1"] <- "Endothelial cell"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "2"] <- "Macrophage"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "3"] <- "B cell"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "4"] <- "Neutrophils"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "5"] <- "dendritic cell"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "6"] <- "NK cell"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "7"] <- "cancer cell1"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "8"] <- "cancer cell2"
combined.obj$`new_cluster`[combined.obj$`new_cluster` == "9"] <- "Plasmacytoid dendritic cell"

combined.obj <- SetIdent(combined.obj, value = combined.obj@meta.data$`new_cluster`)
Seurat::DimPlot(combined.obj, reduction="umap", label = T)  

### 7.visualization ------------------------------------------------------------
features <- c("Cd3e","Ctla4")
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(combined.obj, features = features, ncol = 1)
# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(combined.obj, features = features, ncol = 2)
VlnPlot(combined.obj, features = features, split.by = "orig.ident")
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(combined.obj, features = features)
FeaturePlot(combined.obj, features = features, reduction = "umap",blend = TRUE)
FeaturePlot(combined.obj, features = features, split.by = "orig.ident")
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(combined.obj, features = features) + RotatedAxis()
# Single cell heatmap of feature expression
DoHeatmap(subset(combined.obj, downsample = 100), features = top10$gene, size = 3) + NoLegend()
# Identify conserved cell type markers
sg <- list(bcell = c("Cd19","Cd22","Cd79a","Cd79b"),
           cancer = c("Epcam","Krt18","Krt8","Krt19"),
           dendritic = c("Cd86","Cd83","Clec10a","Ccl22","Cxcl16","Xcr1","Rab7b","Naaa","Btla"),
           endo = c("Pecam1","Tek","Plvap","Thbd","Tmem100","Adgrf5","Podxl","Nostrin","Acvrl1"),
           macrophage = c("Adgre1","Itgam","Itgal","Cd14","Ccl9","F13a1","Cx3cr1","Cd68"),
           neutro= c("Csf3r","Cxcr2S100a8"),
           nkcell = c("Klrd1","Nkg7","Prf1","Gzma","Gzmb","Ccl5","Cma1","Klre1","Klra4"),
           plasma_dend = c("Ccr9","Bst2"),
           tcell = c("Cd3d","Cd5","Cd3e"))
DotPlot(combined.obj, features = sg, cols = c("blue", "red"), dot.scale = 4, assay = "RNA") +
  RotatedAxis()

### 8.計算細胞種類數量 ---------------------------------------------------------
cell_proportions <- combined.obj@meta.data[,c(1,7)]
cell_ratio <- table(cell_proportions) %>% t() %>% as.data.frame()
k_total_cell <- cell_ratio$Freq[1:11] %>% sum()
kl_total_cell <- cell_ratio$Freq[12:22] %>% sum()
k_table <- cell_ratio[1:10,]
kl_table <- cell_ratio[11:22,]
k_table$ratio <- k_table$Freq/k_total_cell
kl_table$ratio <- kl_table$Freq/kl_total_cell
final <- rbind(k_table,kl_table) %>% set_names(c("cluster","ident","freq","ratio"))

# 創建 bar plot
library(ggplot2)
p <- ggplot(data = final, aes(x = ratio, y = ident, fill = cluster)) +
  geom_bar(stat = "identity") + theme_void()
p

# 創建T cell exhausted violin plot
VlnPlot(combined.obj, features = features, ncol = 2)
exh_T <- combined.obj %>%
  subset(., subset = new_cluster == "T cell exhausted")
exh_T <- SetIdent(exh_T, value = exh_T@meta.data$orig.ident) 
VlnPlot(exh_T, features = "Ctla4")

### 9.抓出T cell進行分析 -------------------------------------------------------
act_T <- combined.obj %>%
  subset(., subset = new_cluster == "T cell activated")
act_T <- SetIdent(act_T, value = act_T@meta.data$orig.ident) 
t_de_gene <- FindAllMarkers(act_T) %>% filter(cluster == "kl")
# 創建T cell activated violin plot
VlnPlot(act_T, features = "Tcf7")

# GO分析
library(clusterProfiler)
library(enrichplot)

head(t_de_gene)
organism <- "org.Mm.eg.db"
original_t_gene_list <- t_de_gene$avg_log2FC
names(original_t_gene_list) <- t_de_gene$gene
gsea_t_gene_list <- na.omit(original_t_gene_list) %>% 
  sort(., decreasing = TRUE)

t_gse <- clusterProfiler::gseGO(geneList = gsea_t_gene_list, 
                              ont = "ALL",      # "BP", "MF", "CC"
                              keyType = "SYMBOL", 
                              nPerm = 10000, 
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              pvalueCutoff = 0.05, 
                              verbose = TRUE, 
                              OrgDb = organism, 
                              pAdjustMethod = "none")
# gsea dotplot
dotplot(t_gse, split=".sign") + facet_grid(~.sign)

### 10.抓出cancer cell進行分析 ------------------------------------------------------- 
cancer <- combined.obj %>%
  subset(., subset = new_cluster == c("cancer cell1","cancer cell2"))
cancer <- SetIdent(cancer, value = cancer@meta.data$orig.ident) 
c_de_gene <- FindAllMarkers(cancer) %>% filter(cluster == "kl")

# GO分析
organism <- "org.Mm.eg.db"
original_c_gene_list <- c_de_gene$avg_log2FC
names(original_c_gene_list) <- c_de_gene$gene
gsea_c_gene_list <- na.omit(original_c_gene_list) %>% 
  sort(., decreasing = TRUE)
c_gse <- clusterProfiler::gseGO(geneList = gsea_c_gene_list, 
                              ont = "ALL",      # "BP", "MF", "CC"
                              keyType = "SYMBOL", 
                              nPerm = 10000, 
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              pvalueCutoff = 0.05, 
                              verbose = TRUE, 
                              OrgDb = organism, 
                              pAdjustMethod = "none")
# gsea dotplot
dotplot(c_gse, split=".sign") + facet_grid(~.sign)
