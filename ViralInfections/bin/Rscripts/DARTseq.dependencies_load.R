#############################################################################
# DARTseq.dependencies_load.R     P. Burnham ; 03/26/2018
#############################################################################
## Description
# Set paths an instantiate libraries

#############################################################################
## Dependencies

# Libraries
library(stringr); library(ggplot2); library(data.table); library(caTools);
library(knitr); library(markdown);library(Seurat);library(dplyr);
library(Matrix);library(gdata);library(cowplot); library(scales);

# Working Directory
setwd("/workdir/psb84/DART-seq/ViralInfections/")

# Paths
path = "/workdir/psb84/DART-seq/ViralInfections/"
path.ref = paste0(path, "references/")
path.ref.reo = paste0(path.ref, "Reo/")
path.results = paste0(path, "results/")
path.figs.wdir = paste0(path,"figures/")
path.figs.home = paste0(path,"figures/")
path.stats.home = paste0(path,"stats/")
path.gen.analysis = paste0(path,"general_analysis/")

# Other
cbpal =  c("black","#E69F00", "#0072B2","red")  

#############################################################################
## Functions (alphabetical by function name)

######### Virus Characteristics

# get viral depth by cell and select out for cells passing Seurat QC
virusdepth_qc_cell = function(sample, cell.viability = T){
  vir.depth = data.frame(fread(paste0(path.results,sample,"/",sample,".virdepth.tab"),sep="\t",header = F))
  colnames(vir.depth) = c("Sample","Cell","Segment","Position","Depth")
  if(cell.viability ==T){
    viable.cells = data.frame(fread(paste0(path.results,sample,"/",sample,".viablecells.tab"),sep="\t",header = T))
    vir.depth = vir.depth[vir.depth$Cell %in% as.character(viable.cells$Cell),]
  }
  return(vir.depth)
}

# get viral SNP by cell and select out for cells passing Seurat QC
virusSNP_qc_cell = function(sample, cell.viability = T){
  vir.snp = data.frame(fread(paste0(path.results,sample,"/",sample,".virvar.tab"),sep="\t",header = F))
  colnames(vir.snp) = c("Cell","Segment","Position","Reference","A","C","G","U")
  if(cell.viability == T){
    viable.cells = data.frame(fread(paste0(path.results,sample,"/",sample,".viablecells.tab"),sep="\t",header = T))
    vir.snp = vir.snp[vir.snp$Cell %in% as.character(viable.cells$Cell),]
  }
  return(vir.snp)
}

# aggregate viral depth and produce normalization factors
virus_depth_aggregate = function(sample,cell.viability = T){
  df.bycell = virusdepth_qc_cell(sample)
  df.agg = aggregate(Depth~(Sample+Segment+Position),df.bycell,FUN=sum)
  
  # load in host information
  viable.cells = data.frame(fread(paste0(path.results,sample,"/",sample,".viablecells.tab"),sep="\t",header = T))
  cell.count = nrow(viable.cells)
  umi.count = sum(viable.cells$nUMI)
  
  df.agg$Depth.Norm.Cell = df.agg$Depth / cell.count ;
  df.agg$Depth.Norm.UMI = df.agg$Depth / umi.count ;
  df.agg$Depth.Norm.UMIperCell = df.agg$Depth / (umi.count/cell.count) ;
  
  return(df.agg)  
}

# get major allele fraction from .virvar.tab
determine_majorallele = function(SNP.dataframe){
  nucs = c("A","C","G","U")
  
  agg.A = aggregate(A~(Segment+Position+Reference),SNP.dataframe, FUN = sum)
  agg.C = aggregate(C~(Segment+Position+Reference),SNP.dataframe, FUN = sum)
  agg.G = aggregate(G~(Segment+Position+Reference),SNP.dataframe, FUN = sum)
  agg.U = aggregate(U~(Segment+Position+Reference),SNP.dataframe, FUN = sum)
  
  all.sample = Reduce(function(x, y) merge(x, y,all=TRUE), list(agg.A, agg.C, agg.G, agg.U))
  all.sample$TotalDepth = rowSums(all.sample[4:7])    
  all.sample[all.sample$Reference == "T",]$Reference = "U"
  
  all.sample$Major = nucs[apply(all.sample[,4:7], 1, which.max)]
  all.sample$MajisRef = (all.sample$Major == all.sample$Reference)
  
  all.sample$MAF = sapply(1:nrow(all.sample),
                          function(i) all.sample[i,][[all.sample[i,]$Major]]/(all.sample[i,]$TotalDepth))
  return(all.sample)
}
# get average depth within a region
call.regional.depth = function(df, segment, position, region_size, stat = "Depth" ){
  temp.df = df[df$Segment == segment,]
  temp.depth = temp.df[(temp.df$Position < position) & 
                         (temp.df$Position >=(position - region_size)), ][[stat]]
  return(sum(temp.depth)/region_size)
}

# get average depth within a region
call.regional.depth.trace = function(df, segment, position, region_size, stat = "Depth" ){
  temp.df = df[df$Segment == segment,]
  temp.depth = temp.df[(temp.df$Position < position) & 
                         (temp.df$Position >=(position - region_size)), ][[stat]]
  temp.position = temp.df[(temp.df$Position < position) & 
                            (temp.df$Position >=(position - region_size)), ][["Position"]]
  
  fin.df = merge(data.frame("Position" = ((position-region_size):(position-1)),"Order" = c(region_size:1)), 
                 data.frame("Position" = temp.position, "Stat" = temp.depth),by="Position",all.x=T)
  fin.df$Segment = segment
  if(length(fin.df[is.na(fin.df$Stat),]$Stat) > 0 ){
    fin.df[is.na(fin.df$Stat),]$Stat = 0
  }
  return(fin.df)
}


# produces enrichment by toeholds
toehold_enrichment = function(sample.comp, sample.ref, region_size=200, stat = "Depth.Norm.UMI"){
  probe.stats.short = toehold.props[toehold.props$Sample == sample.comp,]
  
  comp.depth = virus_depth_aggregate(sample.comp)
  ref.depth = virus_depth_aggregate(sample.ref)
  
  probe.stats.short$comp.means = unlist(lapply(c(1:nrow(probe.stats.short)), 
                                               function(x) call.regional.depth(df = comp.depth, 
                                                                               segment = probe.stats.short$Segment[x],
                                                                               position = probe.stats.short$Start[x],
                                                                               region_size = region_size,stat = stat )))
  
  probe.stats.short$ref.means =  unlist(lapply(c(1:nrow(probe.stats.short)), 
                                               function(x) call.regional.depth(df = ref.depth, 
                                                                               segment = probe.stats.short$Segment[x],
                                                                               position = probe.stats.short$Start[x],
                                                                               region_size = region_size,stat = stat )))
  probe.stats.short$stat.used = stat
  probe.stats.short$region.size = region_size
  probe.stats.short$enrichment = probe.stats.short$comp.means / (probe.stats.short$ref.means+10**-10) ;
  return(probe.stats.short)
}

# produces traces by toeholds
toehold_enrichment_trace = function(sample, region_size=200, stat = "Depth.Norm.UMI"){
  if (sample != "Lcell_T3Dpos_Dropseq_0"){
    probe.stats.short = toehold.props[toehold.props$Sample == sample,]
  }else{
    probe.stats.short = toehold.props[toehold.props$Sample == "Lcell_T3Dpos_DARTseq_0",]
  }
  
  sam.depth = virus_depth_aggregate(sample)
  all.trace = call.regional.depth.trace(df =  sam.depth,
                                        segment = probe.stats.short$Segment[1],
                                        position = probe.stats.short$Start[1],
                                        region_size = region_size,stat = stat)
  
  for (i in 2:nrow(probe.stats.short)){
    all.trace = rbind(all.trace,call.regional.depth.trace(df =  sam.depth,
                                                          segment = probe.stats.short$Segment[i],
                                                          position = probe.stats.short$Start[i],
                                                          region_size = region_size,stat = stat))
  }
  
  all.trace$Sample = sample
  all.trace$Measure = stat
  return(all.trace)
}


########## Plotting

plot_rects = function(dfrow,n, xmax,range){
  steps = seq(from = (xmax-range), to = xmax, length.out = (n+1))
  alpha_steps = seq(from =  0, to= 1, length.out = n)
  rect_grad = data.frame(xmin = steps[-(n+1)],
                         xmax = steps[-1],
                         alpha = alpha_steps)
  rect_total = merge(dfrow, rect_grad)
  return(rect_total)
}
