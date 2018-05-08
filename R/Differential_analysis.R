########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

library(Seurat)
library(dplyr)
source("./R/Seurat_functions.R")

#3.1  Compare DE across all major cell types==================
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
#detect changes in gene expression between young and aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in young compared to aged or viceversa. 

# 3.1.1 load data
# Rename ident
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
table(mouse_eyes@ident)
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
new.cluster.ids <- c("0) Pericytes",
                     "1) Retinal pigmented epithelium",
                     "2) Endothelial cells",
                     "3) Pericytes",
                     "4) Endothelial cells",
                     "5) Retinal pigmented epithelium",
                     "6) Smooth muscle cells",
                     "7) Pericytes",
                     "8) Monocytes",
                     "9) Myelinating schwann cells",
                     "10) Endothelial cells",
                     "11) T cells",
                     "12) Monocytes")
mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids,
                                    to = new.cluster.ids)
# 3.1.2 FindAllMarkers.UMI
mouse_eyes_Split <- SplitCells(object = mouse_eyes, split.by = "conditions")
mouse_eyes_aged <- mouse_eyes_Split[[1]]
mouse_eyes_aged.gde <- FindAllMarkers.UMI(object = mouse_eyes_aged)
write.csv(x= mouse_eyes_aged.gde, file="./output/129_B6_aged.csv")

# 3.2 Compare differential expression between subcluster within all major cell types
# plus visualize all major cell types.
#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
# keep the original ident name intact
print("3.2 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# 3.2.1 SubsetData and further split RPE ===============
RPE <- SubsetData(object = mouse_eyes_aged,random.seed = 1,
                  ident.use = new.cluster.ids[grepl("Retinal pigmented epithelium",new.cluster.ids)])
                      
                      
#deselect.cells <- TSNEPlot(object = RPE, do.identify = T)
set.seed(42)
RPE <- FindVariableGenes(object = RPE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
RPE_hv.genes <- head(rownames(RPE@hvg.info), 1000)
RPE <- RunPCA(object = RPE, pc.genes = RPE_hv.genes, pcs.compute = 20, 
              do.print = F)
#PCElbowPlot(object = RPE)
RPE <- FindClusters(object = RPE, reduction.type = "pca", dims.use = 1:5,
                    force.recalc = T, resolution = 0.15, save.SNN = TRUE)
RPE <- RunTSNE(object = RPE, reduction.use = "pca", dims.use = 1:5, 
               do.fast = TRUE)
TSNEPlot(object = RPE, no.legend = TRUE, do.label = TRUE, pt.size = 2,
         label.size = 7)
# generate csv file
RPE.gde <- FindAllMarkers.UMI(object = RPE,test.use = "bimod")
write.csv(x= RPE.gde, file="./output/129_B6_aged_RPE.csv")

# 3.2.2 SubsetData and further split Hematopoietic cells ===============
Hema <- SubsetData(object = mouse_eyes_aged,
                   ident.use = new.cluster.ids[grepl("Monocytes",new.cluster.ids) |
                                               grepl("T cells",new.cluster.ids)])
Hema <- FindVariableGenes(object = Hema, mean.function = ExpMean, dispersion.function = LogVMR, 
                          do.plot = FALSE)
Hema_hv.genes <- head(rownames(Hema@hvg.info), 1000)
set.seed(42)
Hema <- RunPCA(object = Hema, pc.genes = Hema_hv.genes, pcs.compute = 20, do.print = F, 
               pcs.print = 1:5, genes.print = 5)
#PCElbowPlot(object = Hema)
Hema <- FindClusters(object = Hema, reduction.type = "pca", dims.use = 1:5,
                     force.recalc = T,resolution = 1.3, save.SNN = TRUE)
Hema <- RunTSNE(object = Hema, reduction.use = "pca", dims.use = 1:5, 
                do.fast = TRUE)
TSNEPlot(object = Hema, no.legend = TRUE, do.label = TRUE, pt.size = 2,
         label.size = 7)
# generate csv file
Hema.gde <- FindAllMarkers.UMI(object = Hema,test.use = "bimod")
write.csv(x= Hema.gde, file="./output/129_B6_aged_Hematopoietic.csv")

#3.4  Compare DE across all major cell types
#We would need the data for all clusters, as well the subclusters (RPE and hematopoietic cells).
#detect changes in gene expression between 129_B6 and 129_B6_aged, 
#in the different cell types and subtypes. 
#It will also be interesting to check if there is some subtype enriched in 129_B6 compared to 129_B6_aged or viceversa. 
print("3.4 Compare DE across all major cell types")
mouse_eyes.markers <- SplitFindAllMarkers(object = mouse_eyes)

mouse_eyes.gde <- FindAllMarkersbyAge(object = mouse_eyes)
write.csv(x= mouse_eyes.gde, file="./output/mouse_eyes_young_vs_aged.csv")

# 3.5 Compare differential expression between subcluster within all major cell types
# plus visualize all major cell types.
#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
lnames = load(file = "./data/mouse_eyes_alignment.Rda")
lnames
# keep the original ident name intact
print("3.5 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# SubsetData===============
T.cells <- SubsetData(object = mouse_eyes,
                      ident.use = old.ident.ids[(new.cluster.ids %in% "T cells")])
Monocytes <- SubsetData(object = mouse_eyes,
                        ident.use = c(old.ident.ids[(new.cluster.ids %in% "Monocytes")]))
EC <- SubsetData(object = mouse_eyes,
                 ident.use = old.ident.ids[(new.cluster.ids %in% "Endothelial cells")])
RPE <- SubsetData(object = mouse_eyes,
                  ident.use = old.ident.ids[(new.cluster.ids %in% "Retinal Pigment Epithelium")])

Pericytes <- SubsetData(object = mouse_eyes,
                        ident.use = old.ident.ids[(new.cluster.ids %in% "Mesenchymal cells") | 
                                                      (new.cluster.ids %in% "Pericytes") |
                                                      (new.cluster.ids %in% "Smooth muscle cells")])
# split TSNE plot=============
SplitTSNEPlot(T.cells)
SplitTSNEPlot(Monocytes)
SplitTSNEPlot(EC)
SplitTSNEPlot(RPE)
SplitTSNEPlot(Pericytes)


# split and find all markers
table(T.cells@ident)
T.cells.markers <- SplitFindAllMarkers(T.cells, write.csv = TRUE)

table(Monocytes@ident)
Monocytes.markers <- SplitFindAllMarkers(Monocytes, write.csv = TRUE)

table(EC@ident)
EC.markers <- SplitFindAllMarkers(EC, write.csv = TRUE)

table(RPE@ident)
RPE.markers <- SplitFindAllMarkers(RPE, write.csv = TRUE)

table(Pericytes@ident)
Pericytes.markers <- SplitFindAllMarkers(Pericytes, write.csv = TRUE)

#3.6 Compare subcluster between aged vs young===============
print("3.6 Compare subcluster between aged vs young")
# keep the original ident name intact
# Myeloid Cells===============

Myeloid.gde <- FindAllMarkersbyAge(object = Monocytes)
write.csv(x= Myeloid.gde, file="./output/Monocytes_young_vs_aged.csv")

# Perictyes==========
Pericytes.gde <- FindAllMarkersbyAge(object = Pericytes)
write.csv(x= Pericytes.gde, file="./output/Pericytes_young_vs_aged.csv")

# Endothelial Cells==========
table(EC@ident)
Endothelial.gde <- FindAllMarkersbyAge(EC)
write.csv(x= Endothelial.gde, file="./output/Endothelium_young_vs_aged.csv")

# RPE cells=========
RPE.gde <- FindAllMarkersbyAge(RPE)
write.csv(x= RPE.gde, file="./output/RPE_young_vs_aged.csv")
