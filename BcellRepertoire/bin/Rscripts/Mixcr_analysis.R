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
args = commandArgs(trailingOnly = TRUE)
sam=as.character(args[1]) #"9203_7858_68535_HHTNJBGX5_V1-1-500_TTCTGCCT"

## fastq with cell and UMI barcodes
fastq = data.frame(fread(paste0("results/",sam,"/",sam,".tagged.Bcell.reads"),sep = ":",header = F))
fastq2 = data.frame("SeqNum" = c(1:nrow(fastq)),"Cell" = fastq$V10, "UMI" = fastq$V13, "fullID" = paste0(fastq$V10,fastq$V13))

mixcr.results = data.frame(fread(paste0("results/",sam,"/",sam,".tagged.Bcell.alignments.txt"), header = T, sep  = "\t"))
mixcr.results$SeqNum = mixcr.results$readId + 1
full.results = merge(mixcr.results,fastq2, by = "SeqNum", all.x=T)
full.results[is.na(full.results)] = 0

Bcell.list = as.character(as.matrix(fread(paste0("results/",sam,"/",sam,".Bcell.barcodes.txt"),header = F)))

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
exp.matrix = data.frame(fread(paste0("results/",sam,"/",sam,"_expression_matrix.txt"), header = T, sep  = "\t"))
total.counts = data.frame("Cell" = colnames(exp.matrix)[-1],"mRNA" = colSums(data.matrix(exp.matrix[,-1])))

all.data = merge(composite.df,total.counts, by="Cell", all.x = T)
all.data$IGdensity = all.data$both / all.data$mRNA
all.data$play = sample(c("Neither","LC","HC","Both"),replace = T,size = nrow(all.data))
all.data[all.data$light == 0 & all.data$heavy ==0,]$play = "Neither" ;
all.data[all.data$light > 0 & all.data$heavy == 0,]$play = "LC" ;
all.data[all.data$light == 0 & all.data$heavy > 0,]$play = "HC" ;
if(nrow(all.data[all.data$light > 0 & all.data$heavy > 0,])>0){all.data[all.data$light > 0 & all.data$heavy > 0,]$play = "Both" };

write.table(all.data,paste0("results/",sam,"/",sam,".Bcell.IGstats.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("V"),paste0("results/",sam,"/",sam,".Bcell.Vgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("D"),paste0("results/",sam,"/",sam,".Bcell.Dgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("J"),paste0("results/",sam,"/",sam,".Bcell.Jgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)
write.table(region_gen("C"),paste0("results/",sam,"/",sam,".Bcell.Cgenecount.txt"),sep="\t", quote = F, row.names = F, col.names = T)

# see how large the chains are
vdjc = data.frame(cbind(call_section("V"),call_section("D"),call_section("J"),call_section("C")))
vdjc_tab = data.frame(rownames(as.matrix(table(do.call("paste0",vdjc)))),
                   as.numeric(as.matrix(table(do.call("paste0",vdjc)))))
colnames(vdjc_tab) = c("RegionCom","Count")

write.table(vdjc_tab,paste0("results/",sam,"/",sam,".Bcell.chains.txt"),sep="\t", quote = F, row.names = F, col.names = T)

