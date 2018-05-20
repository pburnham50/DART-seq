### Burnham, De Vlaminck: 2018
# Description: Used to calculate highly mutagenic sights and
# analyze transition and transversion rates in single cells.

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")
save_eps=T
expression_min = 2000

# Functions --------------------------------
call_expression = function(sam){
  Sample.data=fread(paste0(path.results,sam,"/",sam,"_expression_matrix.txt"),
                     sep="\t",header = T,stringsAsFactors = F,showProgress = T) ;
  Sample.data <- as.data.frame(Sample.data) ;
  rownames(Sample.data) <- Sample.data$GENE ;
  Sample.data <- Sample.data[,-1] ;
  Sample.data <- Sample.data[,colSums(Sample.data)>=expression_min] ;
  
}
props_in_cluster = function(data.mat, group){
  all.cells = length(grep(group,data.mat@cell.names)) ; 
  cluster.list = unique(as.character(as.matrix(data.mat@ident))) ;
  proplist = sapply(cluster.list,  
                    function(clnum) length(grep(group,SubsetData(data.mat, 
                                                                 ident.use = clnum, 
                                                                 subset.raw = T)@cell.names))) ;
  proplist = as.numeric(unlist(proplist))/all.cells ;
  return(data.frame("Group" = group, "Cluster" =cluster.list, "PropPop"= proplist)) ;
}

# Execution ------------------------------------------

# set four samples
sam1="Lcell_T3Dneg_DARTseq_0" ; sam2="Lcell_T3Dpos_DARTseq_0" ;
sam3="Lcell_T3Dpos_DARTseq_1" ; sam4="Lcell_T3Dpos_Dropseq_0" ;

#call samples individually
Sample1.data <- call_expression(sam1) ;
Sample2.data <- call_expression(sam2) ;
Sample3.data <- call_expression(sam3) ;
Sample4.data <- call_expression(sam4) ;

#create Seurat objects
Sample1 <- CreateSeuratObject(raw.data = Sample1.data, min.cells = 3, min.genes = 250)
Sample2 <- CreateSeuratObject(raw.data = Sample2.data, min.cells = 3, min.genes = 250)
Sample3 <- CreateSeuratObject(raw.data = Sample3.data, min.cells = 3, min.genes = 250)
Sample4 <- CreateSeuratObject(raw.data = Sample4.data, min.cells = 3, min.genes = 250)

#try treating everything as one Seurat and going about analysis
p1 = MergeSeurat(object1 = Sample1, object2 = Sample2, add.cell.id1 = sam1, add.cell.id2 = sam2)
Sample = AddSamples(object = p1, new.data = Sample3.data, add.cell.id = sam3)
Sample = AddSamples(object = Sample, new.data = Sample4.data, add.cell.id = sam4)

# the rest follows directly through the Seurat pipeline
mitogenes = grep(pattern = "mt-", x= rownames(x = Sample@data), value = T)
percent.mito <- Matrix::colSums(Sample@raw.data[mitogenes, ])/Matrix::colSums(Sample@raw.data)
Sample <- AddMetaData(object = Sample, metadata = percent.mito, col.name = "percent.mito")
Sample <- FilterCells(object = Sample, subset.names = c("nGene", "percent.mito"), 
                      low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.025))
Sample <- NormalizeData(object = Sample, normalization.method = "LogNormalize",scale.factor = 10000)
Sample <- FindVariableGenes(object = Sample, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 5, y.cutoff = 0.5)
Sample <- ScaleData(object = Sample, vars.to.regress = c("nUMI", "percent.mito"))
Sample <- RunPCA(object = Sample, pc.genes = Sample@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
host.cells = Sample
host.cells <- RunTSNE(object = host.cells, reduction.use = "pca", 
                      dims.use = 1:4,check_duplicates=F,do.fast = T)
host.cells <- FindClusters(object = host.cells, reduction.type = "pca", 
                      dims.use = 1:4, save.SNN = TRUE,plot.SNN = T,resolution = 0.3)
clustered.tsne <- TSNEPlot(object = host.cells, do.return = TRUE, pt.size = 0.5)

# get a data frame of the cell proportion in each cluster
Sample.past = host.cells
oddcells = c(WhichCells(object = Sample.past, ident = 0),WhichCells(object = Sample.past, ident = 1))
Sample.past = SetIdent(object = Sample.past, cells.use = oddcells, ident.use = 0)
dfcol = rbind(props_in_cluster(data.mat = Sample.past, group=sam1),
              props_in_cluster(data.mat = Sample.past, group=sam2),
              props_in_cluster(data.mat = Sample.past, group=sam3),
              props_in_cluster(data.mat = Sample.past, group=sam4))
subcol = dfcol[dfcol$Cluster %in% c(2,3),]

# get sample markers and export them
tsndf = data.frame(Sample.past@dr$tsne@cell.embeddings,"ID" = Sample.past@ident)
SampleTSNE.markers <- FindAllMarkers(object=Sample,only.pos=F, min.pct=0.10, thresh.use=0.25)

write.table(x = SampleTSNE.markers,file = paste0(path.stats.home,"fourexp.genemarkers.tab"),
            sep = "\t",quote = F, row.names = F, col.names = T)

# Figures -----------------------------------
if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2I.eps"),
                 width=4.05/2.54, height=2.5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot(subcol,aes(Cluster,PropPop))+
  geom_bar(aes(fill = Group),stat="identity",position = "dodge",col="black")+
  ylab("")+xlab("") +
  scale_fill_manual(values=c("black","red","#0072B2","#E69F00"))+
  scale_y_continuous(labels=percent, limits = c(-.01,max(subcol$PropPop)+0.075),
                     breaks = seq(0, max(subcol$PropPop)+0.1,0.05))+
  scale_x_discrete(labels = c("2"="3", "3"="4"))+
  theme_bw()+theme(legend.position="none")+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_text(family="Helvetica",size = 8),
        axis.text.y=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 6))
if(save_eps){dev.off()}