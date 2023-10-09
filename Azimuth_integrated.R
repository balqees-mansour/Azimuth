# https://github.com/satijalab/azimuth/wiki/Generating-an-Azimuth-Reference
# First step is  Azimuth reference preparation to go in Annotation 
# install packages 
update.packages(oldPkgs = c("withr", "rlang"))
remotes::install_github("hhoeflin/hdf5r")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
install.packages("SeuratObject", version = "4.9.9.9060")
remotes::install_github('satijalab/azimuth', ref = 'master')
install.packages(c( "‘hdf5r’, ‘SeuratDisk’"))
devtools::install_version("RcppAnnoy", version = "0.0.16")
install.packages("Annoy")
library(Annoy)
library(Azimuth)
library(hdf5r)
library(SeuratDisk)

# generating Azimuth reference from seurat object 
az_ref <- seurat_ref

#Step 1: Perform SCTransform normalization
az_ref <- SCTransform(az_ref)

# Store the normalized data in the SCT assay
az_ref[["SCT"]] <- az_ref
ref <- AzimuthReference(
  object = az_ref,
  refUMAP = "umap",
  refDR = "pca",
  metadata = c("Cell_label", "Cluster_ID"),
  dims = 1:50,
  k.param = 31
)
saveRDS(ref, file = "ref_Z.Rds")

#Saving the annoy index
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = "idx.annoy")

#==========================================================================================
# The RunAzimuth function can take a Seurat object as input
DimPlot(ref, group.by ="Cell_label", raster=FALSE,
        label = TRUE, repel = TRUE) + NoLegend() 

mm.query <- RunAzimuth(raw.data, 
                       reference = "/home/balqees/Multiple-Myeloma/azimuth_ref")

View(mm.query@meta.data)
mm.query@meta.data$predicted.Cell_label
mm.query@meta.data$predicted_id
DimPlot(mm.query, group.by = "predicted.Cell_label", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()

# Set Idents
Idents(mm.query)

DimPlot(mm.query)

Idents(mm.query) <- "predicted.Cell_label"

DimPlot(mm.query)
DimPlot(ref_z)

# Azimuth normalizes data before mapping, but does not return the results
# normalize the data here before visualization.
mm.query <- NormalizeData(mm.query)

FeaturePlot(mm.query, features = c ("NARF", "SDC1,", "PTP4A3", "COL1A2"))

# Save RDS
saveRDS(mm.query, "../Desktop/Video_Tutorials/Data/lungquery.rds")

ref_z <- readRDS("/home/balqees/Multiple-Myeloma/ref.Rds")
# prepare the query data 
final_query <- readRDS("/home/balqees/Documents/Annot_seurat/final_Downsample_RCPA.Rds")
##================================================================================================================================================
library(Seurat)
library(Azimuth)
library(tidyverse)

# #plot the reference data
rownames(seurat_ref)

view(seurat_ref@meta.data)

DimPlot(seurat_ref, group.by ="Cell_label", raster=FALSE,
        label = TRUE, repel = TRUE) + NoLegend()

DimPlot(lungref, group.by ="ann_finest_level", raster=FALSE,
        label = TRUE, repel = TRUE) + NoLegend()

raw.data #query data 

view(data_query@meta.data)

# The RunAzimuth function can take a Seurat object as input
mm.query <- RunAzimuth(raw.data, 
                       reference = "/home/balqees/Multiple-Myeloma/reference")

view (mm.query@meta.data)

DimPlot(mm.query, group.by = "predicted.celltype.l1", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()

DimPlot(mm.query, group.by = "predicted.celltype.l2", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()
mm.query@meta.data$predicted.celltype.l1
mm.query@meta.data$predicted.celltype.l2
# Set Idents
Idents(mm.query)

DimPlot(mm.query)

Idents(mm.query) <- "predicted.celltype.l2"

DimPlot(mm.query)
DimPlot(ref_z)

# Azimuth normalizes data before mapping, but does not return the results
# normalize the data here before visualization.
mm.query <- NormalizeData(mm.query)

FeaturePlot(mm.query, features = c("SDC1", "SLAMF7", "PTP4A3", "XBP1"))

# Save RDS
saveRDS(mm.query, "../Desktop/Video_Tutorials/Data/lungquery.rds")

ref_z <- readRDS("/home/balqees/Multiple-Myeloma/ref.Rds")
