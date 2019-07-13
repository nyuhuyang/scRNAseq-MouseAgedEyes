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
write.csv(x= mouse_eyes_aged.gde, file="./output/mouse_eyes_aged_129_B6.csv")

# 3.2 Compare differential expression between subcluster within all major cell types
# plus visualize all major cell types.
#http://satijalab.org/seurat/de_vignette.html#perform-de-analysis-using-alternative-tests
# Compare subclusters within aged and young
# keep the original ident name intact
print("3.2 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# 3.2.1 SubsetData and further split RPE ===============
RPE <- SubsetData(object = mouse_eyes,random.seed = 1,
                  ident.use = new.cluster.ids[grepl("Retinal pigmented epithelium",new.cluster.ids)])
                      
                      
#deselect.cells <- TSNEPlot(object = RPE, do.identify = T)
set.seed(42)
RPE <- FindVariableGenes(object = RPE, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
RPE_hv.genes <- head(rownames(RPE@hvg.info), 1000)
RPE <- RunPCA(object = RPE, pc.genes = RPE_hv.genes, pcs.compute = 20, 
              do.print = F)
PCElbowPlot(object = RPE)
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
Hema <- SubsetData(object = mouse_eyes,
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
#lnames = load(file = "./data/mouse_eyes_alignment.Rda")
#lnames
table(mouse_eyes@ident)
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids <- idents$Var1
# keep the original ident name intact
print("3.5 Compare DE between subcluster within all major cell types, and visualize all major cell types")
# SubsetData===============
T.cells <- SubsetData(object = mouse_eyes,
                      ident.use = new.cluster.ids[grepl("T cells",new.cluster.ids)])

Monocytes <- SubsetData(object = mouse_eyes,
                        ident.use = new.cluster.ids[grepl("Monocytes",new.cluster.ids)])

EC <- SubsetData(object = mouse_eyes,
                 ident.use = new.cluster.ids[grepl("Endothelial cells",new.cluster.ids)])
                 
RPE <- SubsetData(object = mouse_eyes,
                  ident.use = new.cluster.ids[grepl("Retinal pigmented epithelium",new.cluster.ids)])
                  
Pericytes <- SubsetData(object = mouse_eyes,
                        ident.use = new.cluster.ids[grepl("Pericytes",new.cluster.ids)])

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

#3.7 Compare subcluster between aged.vs.young(2) and both young(1)===============
mouse_eyes2 <- mouse_eyes
TSNEPlot(object = mouse_eyes2, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 5)+
    ggtitle("TSNE plot of aged and young 129_B6")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle
remove(mouse_eyes)

mouse_eyes_Split <- SplitCells(object = mouse_eyes2, split.by = "conditions")
mouse_eyes_Split[[3]]
young_129_B6_2 <- mouse_eyes_Split[[2]]
p2 <- SplitTSNEPlot(young_129_B6_2, select.plots=1,return.plots=T)
p2[[1]]+ ggtitle("young_129_B6 from aged_and_young_129_B6 dataset")


lnames = load(file = "../scRNAseq-MouseEyes/data/mouse_eyes_alignment.Rda")
lnames
table(mouse_eyes@ident)
idents <- as.data.frame(table(mouse_eyes@ident))
old.ident.ids1 <- idents$Var1
new.cluster.ids1 <- c("0) Pericytes",
                     "1) Pericytes",
                     "2) Endothelial cells",
                     "3) Smooth muscle cells",
                     "4) Retinal pigmented epithelium",
                     "5) Endothelial cells",
                     "6) Pericytes",
                     "7) Pericytes",
                     "8) Hematopoietic cells",
                     "9) Myelinating schwann cells",
                     "10) Endothelial cells",
                     "11) Myelinating schwann cells",
                     "12) Melanocytes")

mouse_eyes@ident <- plyr::mapvalues(x = mouse_eyes@ident,
                                    from = old.ident.ids1,
                                    to = new.cluster.ids1)
TSNEPlot(object = mouse_eyes, no.legend = TRUE, do.label = TRUE,
         do.return = TRUE, label.size = 5)+
    ggtitle("TSNE plot of young 129_B6 and B6")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle
mouse_eyes_Split <- SplitCells(object = mouse_eyes, split.by = "conditions")
mouse_eyes_Split[[3]]
young_129_B6_1 <- mouse_eyes_Split[[1]]
p1 <- SplitTSNEPlot(mouse_eyes, select.plots=1,return.plots=T)
p1[[1]]+ ggtitle("young_129_B6 from 129_B6/B6 dataset")

# 3.7.2 extract certein cell type from young_129_B6_1
# RPE
RPE_1 <- SubsetData(object = young_129_B6_1,random.seed = 1,
                  ident.use = new.cluster.ids1[grepl("Retinal pigmented epithelium",
                                                     new.cluster.ids1)])
TSNEPlot(RPE_1)
select.cells <- WhichCells(object = RPE_1)
deselect.cells <- which(select.cells == "129_B6_CGCTATCTCAGGTTCA")
select.cells <- select.cells[-deselect.cells]
RPE_1 <- SubsetData(RPE_1,cells.use = select.cells)
Cell_RPE_1 <- WhichCells(RPE_1)
colors_1 <- gg_colors(young_129_B6_1,cells.use = Cell_RPE_1) # Extract ColorHexa from Seurat TSNE plot
RPE_1_p <- SplitTSNEPlot(RPE_1, select.plots=1,return.plots=T,
                         pt.size = 5, label.size=15, colors.use = colors_1)
RPE_1_p[[1]]+ ggtitle("RPE in young 129_B6")

# EC
EC_1 <- SubsetData(object = young_129_B6_1,random.seed = 1,
                 ident.use = new.cluster.ids1[grepl("Endothelial cells",
                                                    new.cluster.ids1)])
Cell_EC_1 <- WhichCells(EC_1)
colors_2 <- gg_colors(young_129_B6_1,cells.use = Cell_EC_1) # Extract ColorHexa from Seurat TSNE plot
EC_1_p <- SplitTSNEPlot(EC_1, select.plots=1,return.plots=T,
                        pt.size = 4, label.size=10, 
                        colors.use = colors_2[c(2,1,3)]) # adjust color sequence
EC_1_p[[1]]+ ggtitle("EC in young 129_B6")
FeaturePlot(object = EC_1,reduction.use = "tsne",
                     features.plot = "Ihh", min.cutoff = NA, 
                     cols.use = c("lightgrey","blue"), pt.size = 5)

# Hema
Hema_1 <- SubsetData(object = young_129_B6_1,random.seed = 1,
                   ident.use = new.cluster.ids1[grepl("Hematopoietic cells",new.cluster.ids1)])
Cell_Hema_1 <- WhichCells(Hema_1)
colors_3 <- gg_colors(young_129_B6_1,cells.use = Cell_Hema_1) # Extract ColorHexa from Seurat TSNE plot
Hema_1_p <- SplitTSNEPlot(Hema_1, select.plots=1,return.plots=T,
                        pt.size = 10, label.size=25, 
                        colors.use = colors_3) # adjust color sequence
Hema_1_p[[1]]+ ggtitle("Hematopoietic cells in young 129_B6")

# 3.7.3 extract certein cell type from young_129_B6_2(2)
# RPE
RPE_2 <- SubsetData(object = young_129_B6_2,random.seed = 1,
                    ident.use = new.cluster.ids[grepl("Retinal pigmented epithelium",
                                                       new.cluster.ids)])
TSNEPlot(RPE_2,do.identify = TRUE)
select.cells <- WhichCells(object = RPE_2)
deselect.cells <- which(select.cells == "young_129_B6_ACTGATGGTTGCCTCT")
select.cells <- select.cells[-deselect.cells]
RPE_2 <- SubsetData(RPE_2,cells.use = select.cells)
Cell_RPE_2 <- WhichCells(RPE_2)
colors_RPE_2 <- gg_colors(young_129_B6_2,cells.use = Cell_RPE_2) # Extract ColorHexa from Seurat TSNE plot
RPE_2_p <- SplitTSNEPlot(RPE_2, select.plots=1,return.plots=T,
                         pt.size = 5, label.size=15, colors.use = colors_RPE_2[2:1])
RPE_2_p[[1]]+ ggtitle("RPE in young 129_B6")

# EC
EC_2 <- SubsetData(object = young_129_B6_2,random.seed = 1,
                   ident.use = new.cluster.ids[grepl("Endothelial cells",
                                                      new.cluster.ids)])
Cell_EC_2 <- WhichCells(EC_2)
colors_EC_2 <- gg_colors(young_129_B6_2,cells.use = Cell_EC_2) # Extract ColorHexa from Seurat TSNE plot
EC_2_p <- SplitTSNEPlot(EC_2, select.plots=1,return.plots=T,
                        pt.size = 4, label.size=10, 
                        colors.use = colors_EC_2[c(2,1,3)]) # adjust color sequence
EC_2_p[[1]]+ ggtitle("EC in young 129_B6")
FeaturePlot(object = EC_2,reduction.use = "tsne",
            features.plot = "Ihh", min.cutoff = NA, 
            cols.use = c("lightgrey","blue"), pt.size = 5)
# Hema
Hema_2 <- SubsetData(object = young_129_B6_2,random.seed = 1,
                     ident.use = new.cluster.ids[grepl("Monocytes",new.cluster.ids) |
                                                     grepl("T cells",new.cluster.ids)])
Cell_Hema_2 <- WhichCells(Hema_2)
colors_Hema_2 <- gg_colors(young_129_B6_2,cells.use = Cell_Hema_2) # Extract ColorHexa from Seurat TSNE plot
Hema_2_p <- SplitTSNEPlot(Hema_2, select.plots=1,return.plots=T,
                          pt.size = 10, label.size=10, 
                          colors.use = colors_Hema_2) # adjust color sequence
Hema_2_p[[1]]+ ggtitle("Hematopoietic cells in young 129_B6")

#3.7.4 use ident from both young, color from aged.vs.young ===============
#RPE from 3.7.2

#deselect.cells <- TSNEPlot(object = RPE, do.identify = T)
set.seed(42)
TSNEPlot(RPE_1)
RPE_1 <- FindVariableGenes(object = RPE_1, mean.function = ExpMean, dispersion.function = LogVMR, 
                         do.plot = FALSE)
RPE_hv.genes <- head(rownames(RPE_1@hvg.info), 1000)
RPE_1 <- RunPCA(object = RPE_1, pc.genes = RPE_hv.genes, pcs.compute = 20, 
              do.print = F,seed.use = NULL)
PCElbowPlot(object = RPE_1)
RPE_1 <- FindClusters(object = RPE_1, reduction.type = "pca", dims.use = 1:5,
                    force.recalc = T, resolution = 0.3, save.SNN = TRUE)
RPE_1 <- RunTSNE(object = RPE_1, reduction.use = "pca", dims.use = 1:5, 
               do.fast = TRUE)
TSNEPlot(object = RPE_1, no.legend = TRUE, do.label = TRUE, pt.size = 4,
         do.return = T, label.size = 7) +
    ggtitle("RPE in young 129_B6, re-cluster")+
    theme(text = element_text(size=25),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle


set.seed(42)
TSNEPlot(RPE_2)
RPE_2 <- RunPCA(object = RPE_2, pc.genes = RPE_hv.genes, pcs.compute = 20, 
                do.print = F,seed.use = NULL)
PCElbowPlot(object = RPE_2)
RPE_2 <- FindClusters(object = RPE_2, reduction.type = "pca", dims.use = 1:5,
                      force.recalc = T, resolution = 0.3, save.SNN = TRUE)
RPE_2 <- RunTSNE(object = RPE_2, reduction.use = "pca", dims.use = 1:5, 
                 do.fast = TRUE)
TSNEPlot(object = RPE_2, no.legend = TRUE, do.label = TRUE, pt.size = 4,
         do.return = T, label.size = 7) +
    ggtitle("RPE in young 129_B6, re-cluster")+
    theme(text = element_text(size=25),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle

# EC
EC_1_mat <- as.data.frame(EC_1@dr$tsne@cell.embeddings[,1:2])
EC_1_mat$Cells <- rownames(EC_1_mat)
EC_2_ident <- data.frame(ident= EC_2@ident)
EC_2_ident$Cells <- sub("young_","",rownames(EC_2_ident))
EC_1_matrix <- dplyr::inner_join(EC_1_mat,EC_2_ident,by = "Cells")
ggplot(EC_1_matrix, aes(x=tSNE_1, y=tSNE_2,colour = factor(ident))) + 
    geom_point(size=5) +scale_colour_manual(values=colors_EC_2[c(2,1,3)])+
    ggtitle("EC in young 129_B6")+
    theme(text = element_text(size=20),     #larger text including legend title							
          plot.title = element_text(hjust = 0.5)) #title in middle
