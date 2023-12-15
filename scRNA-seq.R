# 載入packages ---------------------------
library("dplyr")
library("magrittr")
library("Seurat")

# 1.讀檔 -------------------------------------------------------
# 方法1. 讀取10x CellRanger輸出檔 
setwd("/Users/benson/Desktop/Transcriptome/week8") ## 根據自己的需求調整
matrix.dir <- "week8_data"
barcode.path <- file.path(matrix.dir, "barcodes.tsv.gz")
features.path <- file.path(matrix.dir, "features.tsv.gz")
matrix.path <- file.path(matrix.dir, "matrix.mtx.gz")
mat <- Matrix::readMM(file = matrix.path)
feature.names <-  read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names  <-  read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat)  <-  barcode.names$V1
rownames(mat) <-  feature.names$V2
head(feature.names)
class(mat)
dim(mat)
mat[1:10,1:10]

# 方法2: 使用Seurat::ReadMtx 
mat.v2 <- Seurat::ReadMtx(mtx = matrix.path,
                          cells = barcode.path,
                          features = features.path,
                          unique.features = FALSE)  # 建議用FALSE，TRUE會改動名字
mat.v2[1:10, 1:10]
dim(mat.v2)
mat.v3 <- Seurat::ReadMtx(mtx = matrix.path,
                          cells = barcode.path,
                          features = features.path,
                          unique.features = TRUE)  # Make feature names unique(default TRUE)
mat.v3[1:10, 1:10]
dim(mat.v3)
# mat.v2和mat.v3差別在哪裡？ 會有重複的
dplyr::setdiff(rownames(mat.v3), rownames(mat.v2))  # 列出有哪些差異
# 檢查是否有重複的gene name
IDfreqs <- rownames(mat) %>% table()  # 使用table()函数創建頻率表
IDfreqs <- IDfreqs[IDfreqs > 1]  # 抓出重複的
length(IDfreqs)
IDfreqs
feature.names %>% dplyr::filter(V2 == names(IDfreqs)[2]) # 抓出重複第二個
feature.names %>% dplyr::filter(V2 == "Fam205a4") # 抓出"Fam205a4"
# 2. 將基因名子一樣的資料加總起來 --------------------------------
cat(c("<- matrix:", dim(mat), "\n"))
if (length(IDfreqs) > 0){
  for (i in names(IDfreqs)) {
    index <- which(rownames(mat) == i)
    duplicates <- as.matrix(mat[index, ])
    sums <- apply(duplicates, 2, sum, na.rm=TRUE)  # 1: rows, 2: columns
    
    ### 置換成加總數值
    mat[index[1],] <- sums
    
    ### 將重複的資料移除
    for (j in length(index):2) {
      mat <- mat[-index[j], ]
    }
  }
}
cat(c("<- new_matrix:", dim(mat), "\n"))

# 3. 將資料載入Seurat內 -------------------------------------------
seurat.obj <- Seurat::CreateSeuratObject(counts = mat,  ## cells as columns and features as rows
                                         project = "myProject",
                                         min.cells = 10, ## Include features detected in at least this many cells.
                                         min.features = 3, ## Include cells where at least this many features are detected.
                                         assay = "RNA", ## Name of the assay corresponding to the initial input data.
                                         meta.data = NULL)
class(seurat.obj)
seurat.obj

# 4. 熟悉Seurat資料結構 -------------------------------------------
slotNames(seurat.obj)  # return the individual slots in an object.
names(seurat.obj)  # current active assay
Seurat::DefaultAssay(seurat.obj) ### current active assay
seurat.obj@active.assay
seurat.obj[["RNA"]]
out <- Seurat::GetAssayData(seurat.obj,
                            layer = "data" #Specific information to pull (i.e. counts, data, scale.data, ...)
                            ) %>% dim()
out
# 每顆細胞的註解資料，放在meta.data中
seurat.obj[[]] %>% head()
seurat.obj@meta.data %>% head()
seurat.obj[["nCount_RNA"]] %>% head()
seurat.obj@meta.data$nCount_RNA %>% head()

# 5. 計算mitochondrial genes比例，代表死細胞，MT大寫為人 -----------------
# 新增一個Mito_percent column到seurat.obj
seurat.obj[["Mito_percent"]] <- Seurat::PercentageFeatureSet(seurat.obj, pattern = "^mt-")
# 同義seurat.obj@meta.data$Mito_percent <- PercentageFeatureSet(seurat.obj, pattern = "^mt-")
seurat.obj@meta.data %>% head()
seurat.obj[[]] %>% names()

# 6. Cell QC ---------------------------------------
Seurat::VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "Mito_percent"), ncol = 3)
# 移除不要的細胞，利用subset(可以選取符合條件的細胞)
cat(c("-> dimension:", dim(seurat.obj), "\n"))
min_nFeature <- 1000
max_Mito_percent <- 10
seurat.obj <- seurat.obj %>%
  subset(., subset = nFeature_RNA > min_nFeature) %>%
  subset(., subset = Mito_percent < max_Mito_percent)
cat(c("-> dimension:", dim(seurat.obj), "\n"))
Seurat::VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "Mito_percent"), ncol = 3)

# 7. Normalization ---------------------------------------
# Method for normalization.
# LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# CLR: Applies a centered log ratio transformation
# RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. No log-transformation is applied. For counts per million (CPM) set scale.factor = 1e6
seurat.obj <- Seurat::NormalizeData(seurat.obj,
                                    normalization.method = "LogNormalize")
seurat.obj@assays

# 8. Gene selection: high variable genes ---------------------------------------
seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, 
                                           selection.method = "vst", 
                                           nfeatures = 2000 # Number of features to select as top variable features
)
# 看挑出哪些Genes
Seurat::VariableFeatures(seurat.obj) %>% length()
Seurat::VariableFeatures(seurat.obj)[1:20]
# plot of average expression vs. dispersion
hvfinfo <- Seurat::HVFInfo(seurat.obj)  # Get features information for an Assay object.
head(hvfinfo)
hvf <- Seurat::VariableFeatures(seurat.obj)  # Get Highly Variable Features gene list
head(hvf)

hvfinfo %$%  # data frame %$% log(mean) 相當於 log(data frame$mean) 
  plot(log(mean), log(variance), pch = 16, col = "grey90")
hvfinfo[hvf,] %$% 
  points(log(mean), log(variance), pch = 16, col = "red")

# 9. Scale and center the data(PCA之前要做，避免PCA造成誤差，使用Z score) ------
# Before putting the data into PCA for dimensionality reduction we will scale the genes so that they have a mean of 0 and a variance of 1. This is claimed to make the analysis less biased by expression level in the PCA.
seurat.obj <- Seurat::ScaleData(seurat.obj)

# 10. PCA(第一次降維)---------------------------------------
seurat.obj@reductions
seurat.obj <- Seurat::RunPCA(seurat.obj, 
                             features = VariableFeatures(seurat.obj), # Features to compute PCA 
                             npcs = 50, # total Number of PCs to compute and store 
                             verbose = FALSE)
seurat.obj@reductions
Seurat::ElbowPlot(seurat.obj, ndims = 50)  # 抓轉折點（肉眼/數學ex.finfPC package)
Seurat::DimPlot(seurat.obj, reduction="pca")

# [Optional] find the elbow point
pct <- seurat.obj@reductions[["pca"]]@stdev / sum(seurat.obj@reductions[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]  # 哪個PC可以解釋超過90%的variant
co2 <- sort(which((pct[1:length(pct)-1]-pct[2:length(pct)])>0.1), decreasing=T)[1]+1
pcs <- min(co1, co2)
pcs
plot(cumu)

# 11. tSNE & UMAP(第二次降維) -------------------------------------------
seurat.obj <- Seurat::RunTSNE(seurat.obj, 
                              reduction="pca", ## Which dimensional reduction (e.g. PCA, ICA) to use for the tSNE 
                              dims= 1:20, ## Which dimensions to use as input features
                              verbose = FALSE)
seurat.obj@reductions %>% names()
seurat.obj <- Seurat::RunUMAP(seurat.obj, 
                              reduction = "pca", 
                              dims = 1:20,
                              verbose = FALSE)
seurat.obj@meta.data %>% names()
Seurat::FeaturePlot(seurat.obj, reduction="umap", features = "Mito_percent")  # 選擇你想看的gene
DimPlot(seurat.obj, reduction="umap")
DimPlot(seurat.obj, reduction="tsne")

# 12. cluster identification ---------------------------------------
# Nearest-neighbor graph construction 
seurat.obj <- Seurat::FindNeighbors(seurat.obj, 
                                    k.param = 20, ## Defines k for the k-nearest neighbor algorithm
                                    reduction = "pca", 
                                    dims=1:20,
                                    verbose = FALSE
)

seurat.obj <- Seurat::FindClusters(seurat.obj, 
                                   resolution = 1, ## Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
                                   verbose = FALSE
)

seurat.obj[[]] %>% head()
Seurat::Idents(seurat.obj)[1:10]
Seurat::Idents(seurat.obj) <- "RNA_snn_res.1"
Seurat::Idents(seurat.obj)[1:10]
Seurat::DimPlot(seurat.obj, reduction="tsne", group.by = "RNA_snn_res.1")

# 手動更改cluster
seurat.obj@meta.data$new_column <- seurat.obj@meta.data$RNA_snn_res.0.5
seurat.obj@meta.data$new_column[seurat.obj@meta.data$new_column == 5] <- 0
table(seurat.obj@meta.data$new_column)

# 13. Cell doublet prediction -------------------------------------------
# 載入scDblFinder
library("scDblFinder")
count.df <- Seurat::GetAssayData(seurat.obj, slot = "counts", assay = "RNA")
dbl <- scDblFinder(count.df, 
                   clusters = Seurat::Idents(seurat.obj),
                   verbose = FALSE)

seurat.obj[["dbl.class"]] <- dbl$scDblFinder.class
seurat.obj[["dbl.score"]] <- dbl$scDblFinder.score
head(dbl)
seurat.obj[[]] %>% head()
DimPlot(seurat.obj, reduction="tsne", group.by = "dbl.class")

# 14. Cell cycle phase prediction and regression ----------------------------
cc.genes  # 內建cell cycle gene list
# Convert into mouse symbol(非正規)
stringr::str_to_title(cc.genes$s.genes)
mouse.cc.genes <- list(
  s.genes = stringr::str_to_title(cc.genes$s.genes),  # 將字串首字轉為大寫，其他為小寫
  g2m.genes = stringr::str_to_title(cc.genes$g2m.genes)
)
mouse.cc.genes
#### Predict cell cycle phase
seurat.obj <- Seurat::CellCycleScoring(seurat.obj, 
                                       s.features = mouse.cc.genes$s.genes, 
                                       g2m.features = mouse.cc.genes$g2m.genes, 
                                       set.ident = FALSE)
# 觀察是否需要Regression
Seurat::DimPlot(seurat.obj, reduction="umap", group.by = "Phase")

# Regression out of cell cycle phase
regress.obj <- Seurat::ScaleData(seurat.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat.obj)) %>%
  Seurat::RunPCA(npcs = 20, verbose = FALSE) %>%
  Seurat::RunUMAP(reduction = "pca", dims = 1:20, verbose = FALSE) %>%
  Seurat::RunTSNE(reduction="pca", dims= 1:20, verbose = FALSE) %>%
  Seurat::FindNeighbors(reduction = "pca", dims = 1:pcs, verbose = FALSE) %>%
  Seurat::FindClusters(resolution = 0.5, verbose = FALSE)
Seurat::DimPlot(regress.obj, reduction="umap", group.by = "Phase")
