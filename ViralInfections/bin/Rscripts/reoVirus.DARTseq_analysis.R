# reoVirus.DARTseq_analysis.R

#########################################################################################################
# Settings
source("bin/Rscripts/DARTseq.dependencies_load.R")

#####
#Input
#####
args = commandArgs(trailingOnly = TRUE)
sam = as.character(args[1]) #"IS2.Reolow"
expression_min = as.numeric(args[2])

path.sam = paste0(path.results, sam,"/")

transcript.mat= as.data.frame(fread(paste0(path.sam,sam,"_expression_matrix.txt"),
                   sep="\t",header = T,stringsAsFactors = F,showProgress = T))

virus.barcodes = as.data.frame(fread(paste0(path.sam,sam,".totalvircounts.tab"),
                                 sep=" ",header = F,stringsAsFactors = F,showProgress = T))
colnames(virus.barcodes) = c("Vcount","Cell")


#########################################################################################################
# Execution

# Format data
rownames(transcript.mat) <- transcript.mat$GENE
transcript.mat <- transcript.mat[,-1]
transcript.mat <- transcript.mat[,colSums(transcript.mat,na.rm = T)>=expression_min]

# Load matrix into seurat object
Sample <- CreateSeuratObject(raw.data = transcript.mat, min.cells = 3, min.genes = 250)

#remove mitochondrial genes
mitogenes = grep(pattern = "mt-", x= rownames(x = Sample@data), value = T)
percent.mito <- Matrix::colSums(Sample@raw.data[mitogenes, ])/Matrix::colSums(Sample@raw.data)
Sample <- AddMetaData(object = Sample, metadata = percent.mito, col.name = "percent.mito")
Sample <- FilterCells(object = Sample, subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(200, -Inf), high.thresholds = c(Inf, 0.05))

# normalize and find variable genes
Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize",scale.factor = 10000)

Sample <- FindVariableGenes(object = Sample, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 10, y.cutoff = 0.5)

Sample <- ScaleData(object = Sample, vars.to.regress = c("nUMI", "percent.mito"))
Sample <- RunPCA(object = Sample, pc.genes = Sample@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# print cells that make the cut
host.cells = Sample

# find grouping clusters
host.cells <- RunTSNE(object = host.cells, reduction.use = "pca", 
                      dims.use = 1:4,check_duplicates=F,do.fast = T)
host.cells <- FindClusters(object = host.cells, reduction.type = "pca", 
                           dims.use = 1:4, save.SNN = TRUE,resolution = 0.3)
clustered.tsne <- TSNEPlot(object = host.cells, do.return = TRUE, pt.size = 0.5)
SampleTSNE.markers <- FindAllMarkers(object=host.cells,only.pos=F, min.pct=0.10, thresh.use=0.25)


#arange for export
viable_cells = data.frame("Cell" = row.names(host.cells@dr$tsne@cell.embeddings), 
                          "tSNE1" = host.cells@dr$tsne@cell.embeddings[,1],
                          "tSNE2" = host.cells@dr$tsne@cell.embeddings[,2],
                          "Min.Expr.Lvl" = expression_min,
                          "ClusterID" = host.cells@ident,
                          "Sample" = sam)

UMI.mat = data.frame("Cell" =  colnames(host.cells@raw.data),
                           "nUMI" = colSums(host.cells@raw.data))

viable_cells = merge(viable_cells, UMI.mat, by = "Cell",all.x =T)

#########################################################################################################
# Output
#####

## Tables/Files

# publish cell IDs, tSNE cluster coordinates, total expression, viral count
write.table(x = viable_cells,file = paste0(path.sam,sam,".viablecells.tab"),
            quote = F,col.names = T, row.names = F, sep = "\t")
#  clustermarkers
write.table(x = SampleTSNE.markers,file = paste0(path.sam,sam,".cluster.viable.genemarkers.tab"),
            quote = F,col.names = T, row.names = F, sep = "\t")
