#initialize libraries

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

# sample name here
sam="9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA"
wdir=paste0("/workdir/fw262/V5/",sam,"/")
# set workdir above plz
setwd(wdir)

Sample=fread(paste0(wdir,"/",sam,"_expression_matrix.txt"),sep="\t",
             header = T,stringsAsFactors = F,showProgress = T)

Sample <- as.data.frame(Sample) # force to data frame
rownames(Sample) <- Sample$GENE # make the name of rows GENEs
Sample <- Sample[,-1] # take out last column?

Sample.data <- Matrix(as.matrix(Sample),sparse = T)
#Sample.obj <- new("seurat",raw.data=as.matrix(Sample))

Sample <- CreateSeuratObject(raw.data = Sample, min.cells = 2, min.genes = 1, project = "PMBC_MW")

mito.genes <- grep(pattern = "^MT-", x = rownames(x = Sample@data), value = TRUE)
percent.mito <- colSums(Sample@raw.data[mito.genes, ])/colSums(Sample@raw.data)

# no mito gene analysis yet
#Sample <- AddMetaData(object = Sample, metadata = percent.mito, col.name = "percent.mito")
# some plots to check the genes
VlnPlot(object = Sample, features.plot = c("nGene", "nUMI"), nCol = 2)
par(mfrow = c(1, 2))
GenePlot(object = Sample, gene1 = "nUMI", gene2 = "nGene")

Sample <- SubsetData(Sample, subset.name = "nGene", accept.high = 2500)
# normalize the data
Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize", 
                      scale.factor = 10000) # 10 000 is the default given by Seurat

# select for variable genes
Sample <- FindVariableGenes(object = Sample, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality reduction and clustering.
# 
# We can regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a ‘cell-cycle’ score (see example [HERE]) and regress this out as well. In this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content.
# 
# Seurat v2.0 implements this regression as part of the data scaling process. Therefore, the RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData.
Sample <- ScaleData(object = Sample)#, vars.to.regress = "nUMI","percent.mito")

SamplePCA <- RunPCA(object = Sample, pc.genes = Sample@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = SamplePCA, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
PCAPlot(object = SamplePCA, dim.1 = 1, dim.2 = 2)

# In particular PCHeatmap allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets.

PCHeatmap(object = SamplePCA, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

# multiple heat maps here
#PCHeatmap(object = SamplePCA, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

PCElbowPlot(SamplePCA)

################################## create tSNE
SampleTSNE <- FindClusters(object = SamplePCA, reduction.type = "pca", dims.use = 1:10, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

SampleTSNE <- RunTSNE(SampleTSNE, dims.use= 1:10, do.fast = T, check_duplicates = FALSE)
plotTSNEall=TSNEPlot(SampleTSNE, do.return = TRUE, do.label=FALSE)

#################################################################
# analyze tSNE result

# loop through all the clusters and store the cluster expression in a vector
#clusterMarkers<-data.frame(matrix(nrow=nlevels(SampleTSNE@ident),ncol=1)) # initialize data frame variable

clusterMarkers<-vector("list",nlevels(SampleTSNE@ident))
for (clusterNum in (0:(nlevels(SampleTSNE@ident)-1))){
  # the "min.pct" argument requires a gene to be detected at a minimum percentage in either the group in ident.1 or the rest of the cells
  tmp<-FindMarkers(object=SampleTSNE, ident.1=clusterNum, min.pct = 0.25)
  clusterMarkers[clusterNum+1] <- list(tmp) # 
  
}

# find markers for every cluster compared to all remaining cells, report only the positive ones
SampleTSNE.markers <- FindAllMarkers(object=SampleTSNE,only.pos=TRUE, min.pct=0.25, thresh.use=0.25)
# prints a nice list of the result
SampleTSNE.markers %>% group_by(cluster) %>% top_n(2,avg_diff)

# show a violin plot of expression distribution across clusters for Bcell marker
plotClusterBCell=VlnPlot(object = SampleTSNE, features.plot = "MS4A1", do.return = TRUE)
# plot raw UMIs
plotClusterBCellRaw=VlnPlot(object = SampleTSNE, features.plot = "MS4A1", use.raw = TRUE, y.log= TRUE, do.return = TRUE)

# plot expression of MS4A1 marker relative to rest of cells
plotBcellRelativetoRest = FeaturePlot(object=SampleTSNE, features.plot = "MS4A1", do.return = TRUE)

# plot heat map of the top 3 genes for each cluster
top3 <- SampleTSNE.markers %>% group_by(cluster) %>% top_n(3, avg_diff)
plotTop3HeatMap=DoHeatmap(object = SampleTSNE, genes.use = top3$gene, group.by = "ident", remove.key=TRUE)

# plot heat map of the top 3 genes for each cluster
top10 <- SampleTSNE.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
plotTop10HeatMap=DoHeatmap(object = SampleTSNE, genes.use = top10$gene, group.by = "ident", remove.key=TRUE)

# plot_grid(plotTSNEall,plotTSNEall)
# 
# plot1 = TSNEPlot(SampleTSNE, do.return = TRUE) #, no.legend = TRUE, do.label = TRUE)
# plot2 = TSNEPlot(SampleTSNE, do.return = TRUE) #, no.legend = TRUE, do.label = TRUE)
# plot_grid(plot1,plot2)

# show in a 6 panel plot
# par(mfrow=c(2,2))
# plotTSNEall
# plotBcellRelativetoRest
# plotClusterBCell
# plotClusterBCellRaw
# 
# par(mfrow = c(1,2))
# plotTop3HeatMap
# plotTop10HeatMap

#########################################################

# filter for data that in a certain cluster and re-run tSNE clustering to see what's up
########################################
rawDataFilter=SampleTSNE@raw.data[SampleTSNE@ident==3] # number here corresponds to cluster you want to filter for
cell_barcodes <-colnames(rawDataFilter) # contains 12 cell barcode for cells in cluster of interest
setwd(wdir)
write(x=cell_barcodes, file="Bcell_cellbarcodes.txt") ### file called "Bcell_cellbarcodes.txt" contains modifiled cell barcodes
####################################################################

# create Seurat Object from filtered data
SampleFiltered=CreateSeuratObject(raw.data = rawDataFilter, min.cells = 2, min.genes = 1, project = "Filtered Data")

SampleFiltered <- SubsetData(SampleFiltered, subset.name = "nGene", accept.high = 2500)
# normalize the data
SampleFiltered <- NormalizeData(object = SampleFiltered, normalization.method = "LogNormalize", 
                        scale.factor = 10000) # 10 000 is the default given by Seurat

# select for variable genes
SampleFiltered <- FindVariableGenes(object = SampleFiltered, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

SampleFiltered <- ScaleData(object = SampleFiltered, vars.to.regress = "nUMI")
SampleFiltered <- RunPCA(object = SampleFiltered, pc.genes = SampleFiltered@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

SampleFiltered <- FindClusters(object = SampleFiltered, reduction.type = "pca", dims.use = 1:10, 
                           resolution = 0.6, print.output = 0, save.SNN = TRUE)

SampleFiltered <- RunTSNE(SampleFiltered, dims.use= 1:10, do.fast = T, check_duplicates = FALSE)
TSNEPlot(SampleFiltered, do.return = TRUE)

SampleFiltered.markers <- FindAllMarkers(object=SampleFiltered,only.pos=TRUE, min.pct=0.25, thresh.use=0.25)

top10 <- SampleFiltered.markers %>% group_by(cluster) %>% top_n(10, avg_diff)
plotTop10HeatMap=DoHeatmap(object = SampleFiltered, genes.use = top10$gene, group.by = "ident", remove.key=TRUE)


# store clustering result of resolution 0.6
#SampleTSNE <- StashIdent(object=SampleTSNE, save.name="ClusterNames_0.6Res")

#SampleTSNE <- FindClusters(object = SamplePCA, reduction.type = "pca", dims.use = 1:10, resolution = 1.2, print.output = 0, save.SNN = TRUE)

#SampleTSNE <- RunTSNE(SampleTSNE, dims.use= 1:10, do.fast = T, check_duplicates = FALSE)

#plot1 <- TSNEPlot(object = SampleTSNE, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
#plot2 <- TSNEPlot(object = SampleTSNE, do.return = TRUE, group.by = "ClusterNames_0.6Res", no.legend = TRUE, do.label = TRUE)
#plot_grid(plot1, plot2)



# Phil's code here for extracting B cell info
#Bcells = rownames(data.frame(SampleTSNE@ident[SampleTSNE@ident==top10[top10$gene == "MS4A1",]$cluster]))
#nonBcells = rownames(data.frame(SampleTSNE@ident[SampleTSNE@ident!=top10[top10$gene == "MS4A1",]$cluster]))

#write.table(Bcells, paste0("/workdir/psb84/Drop_TC/V1/",sam,"/",sam,".confirmedBcell.list"),quote = F,row.names = F,col.names = F)
#write.table(nonBcells, paste0("/workdir/psb84/Drop_TC/V1/",sam,"/",sam,".nonBcell.list"),quote = F,row.names = F,col.names = F)

allPlot<-TSNEPlot(SampleTSNE, do.return = TRUE, do.label=FALSE,no.legend=TRUE)
FeaturePlot(object=SampleTSNE, features.plot = c("IL7R","CD8A","MS4A1"), cols.use = c("grey", "blue"), reduction.use = "tsne")

