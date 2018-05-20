### Burnham, De Vlaminck: 2018
# Description: Used to calculate subgroups in cell monoculture following
# infection and to identify average viral transcript abundance in each
# single cell

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")

args = commandArgs(trailingOnly = TRUE) 
sam = as.character(args[1]) # Default="Lcell_T3Dpos_DARTseq_0"
min.genes = as.numeric(args[2]) #default=2000, minimum number of genes
save_eps=as.logical(args[3]) 
path.sam = paste0(path.results, sam,"/")

# Execution ---------------------------------------------

# call expression matrix
transcript.mat= as.data.frame(fread(paste0(path.sam,sam,"_expression_matrix.txt"),
                   sep="\t",header = T,stringsAsFactors = F,showProgress = T))
virus.barcodes = as.data.frame(fread(paste0(path.sam,sam,".totalvircounts.tab"),
                                 sep=" ",header = F,stringsAsFactors = F,showProgress = T))
colnames(virus.barcodes) = c("Count","Cell")

# Format data
rownames(transcript.mat) <- transcript.mat$GENE
transcript.mat <- transcript.mat[,-1]
transcript.mat <- transcript.mat[,colSums(transcript.mat,na.rm = T)>=min.genes]

# continue as following Seurat pipeline
Sample <- CreateSeuratObject(raw.data = transcript.mat, min.cells = 3, min.genes = 150)
mitogenes = grep(pattern = "mt-", x= rownames(x = Sample@data), value = T)
percent.mito <- Matrix::colSums(Sample@raw.data[mitogenes, ])/Matrix::colSums(Sample@raw.data)
Sample <- AddMetaData(object = Sample, metadata = percent.mito, col.name = "percent.mito")
Sample <- FilterCells(object = Sample, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.05))
Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize",scale.factor = 10000)
Sample <- FindVariableGenes(object = Sample, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 0.5)
Sample <- ScaleData(object = Sample, vars.to.regress = c("nUMI", "percent.mito"))
Sample <- RunPCA(object = Sample, pc.genes = Sample@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# determinine clusters and differentially expressed genes in each cluster
host.cells = Sample
host.cells <- RunTSNE(object = host.cells, reduction.use = "pca", 
                      dims.use = 1:5,check_duplicates=F,do.fast = T)
host.cells <- FindClusters(object = host.cells, reduction.type = "pca", 
                            dims.use = 1:5, save.SNN = TRUE,plot.SNN = T,resolution = 0.3)
clustered.tsne <- TSNEPlot(object = host.cells, do.return = TRUE, pt.size = 0.5)
SampleTSNE.markers <- FindAllMarkers(object=host.cells,only.pos=F, min.pct=0.10, thresh.use=0.25)

#calculate top 4 genes in each cluster
tops = SampleTSNE.markers %>% group_by(cluster) %>% top_n(4, avg_logFC)

# find markers for every cluster compared to all remaining cells, report only the positive ones
df.collapse = data.frame("Gene" = rownames(host.cells@raw.data),"Exp"=log10(rowSums(host.cells@raw.data)))
SampleTSNE.highpval.markers <- SampleTSNE.markers[SampleTSNE.markers$p_val_adj < 10**-17,]

#Determine cells from above containing virus
tlb = data.frame(host.cells@dr$tsne@cell.embeddings)
tlb$Cell = row.names(tlb)

cluster.identity = data.frame(host.cells@ident)
cluster.identity$Cell = row.names(cluster.identity)
colnames(cluster.identity)[1] = "ClusterID"

raw.rna.data = host.cells@raw.data
total.exp.cell = data.frame("Cell" = colnames(raw.rna.data),
                            "Total.exp" = colSums(raw.rna.data))

vhost = merge(merge(merge(tlb,virus.barcodes,by="Cell",all.x=T),
                    cluster.identity,by="Cell",all.x=T),
                    total.exp.cell,by="Cell",all.x=T)
vhost[is.na(vhost$Count),]$Count = 0
vhost$FracViral = vhost$Count / ( vhost$Total.exp )
vhost$Sample = sam

# average viral fraction per cluster
viral.clusters = aggregate(FracViral ~ ClusterID,vhost, FUN = mean)
colnames(viral.clusters)[1] = "cluster"

# statistical tests
cl =3 ## cluster identified manually
wilcox.test(vhost[vhost$ClusterID == 2,]$FracViral,vhost[(vhost$ClusterID != 2 ),]$FracViral)$p.value

# Tables  -------------------------------------------------------

write.table(x = vhost, file = paste0(path.stats.home,sam,".viruscounts.clusters.tab"),
            sep = "\t", quote = F, row.names = F, col.names = T) ;

# Figures -------------------------------------------------------

if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2G.eps"), 
                 width=5/2.54, height=5.5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=4)}

DoHeatmap(object = host.cells, genes.use = tops$gene, slim.col.label = T, 
          remove.key = T,cex.row = 7,cex.col = 0)
if(save_eps){dev.off()}

if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2H.eps"), 
                 width=4.06/2.54, height=2.4/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=4)}

ggplot(vhost,aes(as.character(ClusterID),FracViral))+
  geom_boxplot(outlier.colour = NA)+
  scale_x_discrete(labels = c("2"="3", "3"="4", "0"="1", "1"="2"))+
  scale_y_continuous(limits = c(-0.0001,0.0031),
                     breaks = seq(0,0.003,0.001),
                     labels = scales::percent)+
  theme_bw()+xlab("")+ylab("Viral transcripts     ")+
  theme(legend.position="none")+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_blank(),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 6))

if(save_eps){dev.off()}