#data link https://www.dropbox.com/s/79q6dttg8yl20zg/immune_alignment_expression_matrices.zip?dl=1

library(Seurat)
library(harmony)
#Read10x can be used if the output matrix comes from 10x cell ranger. here we obtain the matrix file itself as a single file.

ctrl.data <- read.table(file = "C:/Users/immune_alignment_expression_matrices/immune_control_expression_matrix.txt.gz", sep = "\t")
stim.data <- read.table(file = "C:/Users/immune_alignment_expression_matrices/immune_stimulated_expression_matrix.txt.gz", sep = "\t")

ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 3)

ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")

VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)

VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




#########

stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 3)

stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")

VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

stim <- subset(stim, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & nCount_RNA < 6000 & percent.mt < 5)

VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




############################ integrate ###

immune.merged <- merge(ctrl,stim)

immune.merged <- NormalizeData(immune.merged, verbose = FALSE)
immune.merged <- FindVariableFeatures(immune.merged)
immune.merged <- ScaleData(immune.merged, verbose = FALSE)
immune.merged <- RunPCA(immune.merged, npcs = 30, verbose = FALSE)


immune.integrated <- RunHarmony(immune.merged, "orig.ident")

DimPlot(object = immune.integrated, reduction = "harmony", pt.size = .1, group.by = "orig.ident")

immune.integrated <- FindNeighbors(immune.integrated, reduction = "harmony",dims = 1:30)
immune.integrated <- FindClusters(immune.integrated, resolution = 1)


immune.integrated <- RunUMAP(immune.integrated, reduction = "harmony",  dims = 1:30)


DimPlot(immune.integrated, reduction = "umap", group.by = "orig.ident", pt.size = .1) # compare with merged data only
DimPlot(immune.integrated, reduction = "umap", label = TRUE,  pt.size = .1)



####################

DefaultAssay(immune.integrated) <- 'RNA'

## check marker gene expressions

FeaturePlot(immune.integrated, features = c("CD3D", "GNLY", "IFI6"), 
            split.by = "orig.ident", max.cutoff = 3, cols = c("grey","red"), reduction = "umap")


FeaturePlot(immune.integrated, features = "GNLY")

## find markers of clusters
out_path <- "D:/sc_analysis/"
immune.cluster.markers <- FindAllMarkers(immune.integrated,logfc.threshold = 0.1,only.pos = T)

write.csv(immune.cluster.markers,paste(out_path,'immune_cluster_markers.csv',sep=''))
### use marker genes to annotate (clusters to cell types)

#just dummy, use your own annotation
immune.integrated <- RenameIdents(immune.integrated, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")
immune.integrated$CellTypes <- Idents(immune.integrated)

DimPlot(immune.integrated, reduction = "umap", group.by = "CellTypes", pt.size = .1) 


DimPlot(immune.integrated, label = TRUE)

ct_markers <- FindAllMarkers(immune.integrated,logfc.threshold = 0.1,only.pos = T)
