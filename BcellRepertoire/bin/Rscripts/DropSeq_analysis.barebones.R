#####
#initialize libraries

library(Seurat) ; library(dplyr) ; library(Matrix) ; library(data.table)

#####
# Settings

args = commandArgs(trailingOnly = TRUE)
sam=as.character(args[1]) #"9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA"
sam.path = paste0("results/",sam,"/")

plusgene = as.logical(args[2]) #use to add samples that express MS4A1 but may not be in Bcell cluster

#####
# Read in tables

Sample=fread(paste0(sam.path,sam,"_expression_matrix.txt"),sep="\t",header = T,stringsAsFactors = F,showProgress = T)

Sample <- as.data.frame(Sample) # force to data frame
rownames(Sample) <- Sample$GENE # make the name of rows GENEs
Sample <- Sample[,-1] # take out first column

Sample.data <- Matrix(as.matrix(Sample),sparse = T)
Sample <- CreateSeuratObject(raw.data = Sample, min.cells = 2, min.genes = 1, project = "PMBC_MW")
Sample <- SubsetData(Sample, subset.name = "nGene", accept.high = 2500)

# normalize the data
Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize", scale.factor = 10000) # 10 000 is the default given by Seurat

# select for variable genes
Sample <- FindVariableGenes(object = Sample, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

# Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality reduction and clustering.
# We can regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a ‘cell-cycle’ score (see example [HERE]) and regress this out as well. In this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content.
# Seurat v2.0 implements this regression as part of the data scaling process. Therefore, the RegressOut function has been deprecated, and replaced with the vars.to.regress argument in ScaleData.
Sample <- ScaleData(object = Sample)#, vars.to.regress = "nUMI","percent.mito")

SamplePCA <- RunPCA(object = Sample, pc.genes = Sample@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

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

# calculates the cluster corresponding to the Bcell population
BCell_clustnum = as.numeric(as.matrix(SampleTSNE.markers[rownames(SampleTSNE.markers) == "MS4A1",]$cluster))


#########################################################

# filter for data that in a certain cluster and re-run tSNE clustering to see what's up
########################################
rawDataFilter=SampleTSNE@raw.data[SampleTSNE@ident== BCell_clustnum] # number here corresponds to cluster you want to filter for

# contains 12 cell barcode for cells in cluster of interest
cell_barcodes_cluster <- colnames(rawDataFilter) 

if(plusgene == T){
  sls = t(SampleTSNE@raw.data[rownames(SampleTSNE@raw.data) == "MS4A1",] )
  cell_barcodes_gene <- rownames(as.matrix(sls[sls>0,]))
  all_barcodes = intersect(cell_barcodes_cluster, cell_barcodes_gene) #unique(c(cell_barcodes_gene, cell_barcodes_cluster))
}else{all_barcodes = cell_barcodes_cluster}

write(x= all_barcodes, file=paste0(sam.path,sam,".Bcell.barcodes.txt")) ### file called "Bcell_cellbarcodes.txt" contains modifiled cell barcodes


cell.names = as.character(as.matrix(SampleTSNE@cell.names)) #data.frame()
cell.embeddings = SampleTSNE@dr$tsne@cell.embeddings #data.frame()
cell.cluster = SampleTSNE@ident
cell.UMI = data.frame("Cell"=colnames(SampleTSNE@raw.data),"nUMI"=SampleTSNE@meta.data$nUMI)

cell.descriptions = data.frame("Cell" = cell.names,"Cluster" = cell.cluster, 
                               "tSNE_1" = cell.embeddings[,1],"tSNE_2" = cell.embeddings[,2])

final.desc= merge(merge(cell.descriptions, data.frame("Cell" = all_barcodes, "Bcell" = 1), all.x = T),cell.UMI,all.x=T)
final.desc[is.na(final.desc$Bcell),]$Bcell = 0

write.table(x= final.desc, file=paste0(sam.path,sam,".allcells.descriptions.txt"), quote = F, sep = "\t",col.names = T, row.names = F) 
