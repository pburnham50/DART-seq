### Bcell_IG_analysis.R

#####
# Description
#
# Use to take MIXCR results from a data set and,
#   1) compare them to Cell and UMI Barcode information,
#   2) get stats on clonotypes for each cell,
#   3) determine cells with both heavy and light chain,
#   4) overlay on top of tSNE.
#
#####
# Libraries
library(ggplot2); library(data.table); library(caTools);

#####
# Functions

heavy_light_presence = function(CellID = "ACGCCTGTACTG"){
  
  cell.df = df[df$Cell == CellID,]
  
  if(nrow(cell.df)>0){
    cell.df$V3 = substr(cell.df$bestVGene,1,3) ;
    cell.df$D3 = substr(cell.df$bestDGene,1,3) ;
    cell.df$J3 = substr(cell.df$bestJGene,1,3) ;
    cell.df$C3 = substr(cell.df$bestCGene,1,3) ;
    cell.df$IGL = as.numeric((cell.df$V3 == "IGL")|(cell.df$D3 == "IGL")|(cell.df$J3 == "IGL")|(cell.df$C3 == "IGL"))
    cell.df$IGK = as.numeric((cell.df$V3 == "IGK")|(cell.df$D3 == "IGK")|(cell.df$J3 == "IGK")|(cell.df$C3 == "IGK"))
    cell.df$IGH = as.numeric((cell.df$V3 == "IGH")|(cell.df$D3 == "IGH")|(cell.df$J3 == "IGH")|(cell.df$C3 == "IGH"))
    
    cell.frame = data.frame("Cell" = CellID)
    cell.frame$numIGL = sum(cell.df$IGL)
    cell.frame$numIGK = sum(cell.df$IGK)
    cell.frame$numIGH = sum(cell.df$IGH)
  }else{cell.frame = data.frame("Cell" = as.character(CellID), "numIGL" = 0, "numIGK" = 0, "numIGH" = 0 )}
  return(as.matrix(cell.frame))
}

region_gen = function(gene="C"){
  slodf = df[[paste0("best",gene,"Gene")]]
  meddf = data.frame(rownames(as.matrix(table(slodf))),as.numeric(as.matrix(table(slodf))))
  colnames(meddf) = c(paste0(gene,"gene"),paste0(gene,"genecount"))
  return(meddf)
}

call_section = function(gene="C"){
  slodf = df[[paste0("best",gene,"HitScore")]]
  truth = as.character(as.numeric(slodf>0))
  return(gsub("1",gene,gsub("0","",truth)))
}


#####
# Import data
args = c(sample.list[2],F,T)#commandArgs(trailingOnly = TRUE)
#####
sam=as.character(args[1]) #"9203_7858_68535_HHTNJBGX5_V1-1-500_TTCTGCCT"

## fastq with cell and UMI barcodes
fastq = data.frame(fread(paste0("V5/",sam,"/",sam,".tagged.Bcell.reads"),sep = ":",header = F))
fastq2 = data.frame("SeqNum" = c(1:nrow(fastq)),"Cell" = fastq$V10, "UMI" = fastq$V13, "fullID" = paste0(fastq$V10,fastq$V13))

mixcr.results = data.frame(fread(paste0("V5/",sam,"/",sam,".tagged.Bcell.alignments.txt"), header = T, sep  = "\t"))
mixcr.results$SeqNum = mixcr.results$readId + 1
full.results = merge(mixcr.results,fastq2, by = "SeqNum", all.x=T)
full.results[is.na(full.results)] = 0

Bcell.list = as.character(as.matrix(fread(paste0("V5/",sam,"/",sam,".Bcell.barcodes.txt"),header = F)))

#####
# Execution
full.results$totalscore = full.results$bestVHitScore + full.results$bestDHitScore + full.results$bestJHitScore + full.results$bestCHitScore
full.results$allregions = (full.results$bestVGene != "")&(full.results$bestDGene != "")&(full.results$bestJGene != "")&(full.results$bestCGene != "")

##remove duplicates based on best score
df = full.results[order(full.results[,'fullID'],-full.results[,'totalscore']),]
df = df[!duplicated(df$fullID),]

composite.df = data.frame(matrix(unlist(lapply(Bcell.list, function(x) heavy_light_presence(x))),ncol = 4, byrow = T))
colnames(composite.df) = c("Cell","IGL","IGK","IGH")
composite.df$light = as.numeric(as.matrix(composite.df$IGK)) + as.numeric(as.matrix(composite.df$IGL))
composite.df$heavy = as.numeric(as.matrix(composite.df$IGH)) 
composite.df$both = composite.df$light +composite.df$heavy

###### read in expression matrix to give total counts
exp.matrix = data.frame(fread(paste0("V5/",sam,"/",sam,"_expression_matrix.txt"), header = T, sep  = "\t"))
total.counts = data.frame("Cell" = colnames(exp.matrix)[-1],"mRNA" = colSums(data.matrix(exp.matrix[,-1])))

all.data = merge(composite.df,total.counts, by="Cell", all.x = T)
all.data$IGdensity = all.data$both / all.data$mRNA
all.data$play = sample(c("Neither","LC","HC","Both"),replace = T,size = nrow(all.data))
all.data[all.data$light == 0 & all.data$heavy ==0,]$play = "Neither" ;
all.data[all.data$light > 0 & all.data$heavy == 0,]$play = "LC" ;
all.data[all.data$light == 0 & all.data$heavy > 0,]$play = "HC" ;
if(nrow(all.data[all.data$light > 0 & all.data$heavy > 0,])>0){all.data[all.data$light > 0 & all.data$heavy > 0,]$play = "Both" };

write.table(all.data,paste0("V5/",sam,"/",sam,".Bcell.IGstats.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("V"),paste0("V5/",sam,"/",sam,".Bcell.Vgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("D"),paste0("V5/",sam,"/",sam,".Bcell.Dgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("J"),paste0("V5/",sam,"/",sam,".Bcell.Jgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("C"),paste0("V5/",sam,"/",sam,".Bcell.Cgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)

# see how large the chains are
vdjc = data.frame(cbind(call_section("V"),call_section("D"),call_section("J"),call_section("C")))
vdjc_tab = data.frame(rownames(as.matrix(table(do.call("paste0",vdjc)))),
                   as.numeric(as.matrix(table(do.call("paste0",vdjc)))))
colnames(vdjc_tab) = c("RegionCom","Count")

write.table(vdjc_tab,paste0("V5/",sam,"/",sam,".Bcell.chains.txt"),sep="\t", quote = F, row.names = F, col.names = T)

#####################################################
if(as.logical(args[3]) == T){
  
  library(Seurat) ; library(dplyr) ; library(Matrix) ; library(data.table)
  
  #####
  # Settings
  
  main.path = "/workdir/psb84/Drop_BCells/" 
  fig.path = "/home/psb84/"
  sam=as.character(args[1]) #"9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA"
  sam.path = paste0(main.path,"V5/",sam,"/")
  
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
  
  
  
  tsne.df = data.frame("Cell" = rownames(plotTSNEall$data),plotTSNEall$data)
  acno = merge(tsne.df,all.data, by ="Cell", all=T)

  pdf(file=paste0(path.figs,sam,".Bcell.target.highlight.eps"), width=14/2.54, height=9/2.54, paper="special", bg="white",
      fonts="Helvetica",colormode="cmyk",pointsize=1)
  
  ggplot(acno, aes(x = x, y = y)) + 
    geom_point(data = subset(acno[is.na(acno$play),]), alpha = 0.5, col="grey",fill="grey",shape=20) +
    geom_point(data = subset(acno[!is.na(acno$play),]), aes(col = play), alpha = 0.7,shape=20)+
    scale_color_manual(values = c("red","blue","darkgreen","grey"))+
    theme(legend.position="none")+
    xlab("tSNE-1") + ylab("tSNE-2") +
    theme(axis.title.x=element_text(family="Helvetica", size = 8),
          axis.title.y=element_text(family="Helvetica", size = 8),
          axis.text=element_text(family="Helvetica",size = 8),
          plot.title=element_text(family="Helvetica",size = 8))
  dev.off()
  
}
