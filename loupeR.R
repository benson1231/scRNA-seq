library(Seurat)
library(remotes)
library(hdf5r)
library(magrittr)

remotes::install_github("10XGenomics/loupeR")
loupeR::setup()

# import the library
library(loupeR)

setwd("/Users/benson/Desktop/Transcriptome/week12")
seurat_obj <- file.path("test_object.RDS") %>% readRDS()
# convert the SeuratObject named `seurat_obj` to a Loupe file
loupeR::create_loupe_from_seurat(seurat_obj, output_name = "test_object")
