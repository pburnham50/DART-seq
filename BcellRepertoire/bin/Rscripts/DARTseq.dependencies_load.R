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

# Paths
path = "./"
path.ref = paste0(path, "references/")
path.results = paste0(path, "results/")
path.figs.home = paste0(path,"figures/")
path.stats.home = paste0(path,"stats/")

# Other
cbpal =  c("red","#56B4E9","#E69F00","black")  

#############################################################################
## Functions (alphabetical by function name)


cell_frac = function(mrna.bin,dframe, look = "Both"){
  df = dframe ;
  df[df$Binexp == mrna.bin,] ;
  frac = nrow(df[(df$Binexp == mrna.bin) & (df$play == look | df$play == "Both"),])/nrow(df[df$Binexp == mrna.bin,]) ;
  nce = nrow(df[df$Binexp == mrna.bin,]) ;
  if (look == "Either"){
    frac = nrow(df[(df$Binexp == mrna.bin) & (df$play != "Neither"),])/nrow(df[df$Binexp == mrna.bin,]) ;
    nce = nrow(df[df$Binexp == mrna.bin,]) ;
  }
  return(cbind(nce,frac)) ;
}


quick_IG_stats = function(sam,low.RNA = 1200){
  df = data.frame(fread(paste0(path.results,sam,"/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
  df.sub = df[df$mRNA >= low.RNA, ]
  df.hits = data.frame(table(df.sub$play))
  
  either.frac = sum(df.hits[df.hits$Var1 != "Neither",]$Freq)/sum(df.hits$Freq)
  both.frac = sum(df.hits[df.hits$Var1  == "Both",]$Freq)/sum(df.hits$Freq)
  
  stat.tab = cbind(c("Sample","lower UMI", "Fraction Either", "Fraction Both"),
                   c(sam, low.RNA, round(either.frac,4), round(both.frac,4)))
  
  return(stat.tab)
}

binCells = function(sam, binlevel = 250, max.count = 4000, look = "LC", cells.per = 10){
  df = data.frame(fread(paste0(path.results,sam,"/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
  df$Binexp = cut(df$mRNA,seq(0,max.count,binlevel),labels = F)*binlevel ;
  levs = seq(0,max.count,binlevel) ;
  finale = sapply(levs, function(x) cell_frac(mrna.bin = x, df, look = look))
  nce = data.frame(levs, t(finale), sam, look)   
  colnames(nce) = c("BinExp","nCell","propCell", "Sample", "IGtype")
  return(nce[nce$nCell >= cells.per,])
}


sigmoid_fit = function(sample, bin.lvl, cells.per = 10, max.count, look, setb=T,param.b=10){
  
  alc = binCells(sam = sample,binlevel = bin.lvl, max.count = max.count , look = look, cells.per = cells.per ) ;
  x.exp = alc$BinExp ;
  y.exp = alc$propCell ;
  
  # best fit
  if(setb==T){
    y.thr = nls(y.exp~1/(1+exp(-param.b*(x.exp-c))),start = list(c=1500),
                control = nls.control(maxiter = 4000, minFactor = 1/10000,warnOnly = T)) ;
  }else{
    y.thr = nls(y.exp~1/(1+exp(-b*(x.exp-c))),start = list(b=1/1000,c=1500),
                control = nls.control(maxiter = 4000, minFactor = 1/10000,warnOnly = T)) ;
  }
  # get summary
  df.summary = summary(y.thr)
  
  return(data.frame(df.summary$coefficients,"Sample" = sample))
}

sigmoid_fit_produce = function(df, max.value=20000, step.size = 100,bval){
  
  # coefficients of fit
  b.fit=bval ;
  c.fit.mean=df[rownames(df) == "c",]$Estimate
  c.fit.low=df[rownames(df) == "c",]$Estimate - df[rownames(df) == "c",]$Std..Error ;
  c.fit.high=df[rownames(df) == "c",]$Estimate + df[rownames(df) == "c",]$Std..Error ;
  
  xaz = seq(0,max.value,step.size) ;
  yaz.low = 1/(1+exp(-b.fit*(xaz-c.fit.low))) ;
  yaz.mid = 1/(1+exp(-b.fit*(xaz-c.fit.mean))) ;
  yaz.high = 1/(1+exp(-b.fit*(xaz-c.fit.high))) ;
  
  return(data.frame("BinExp" = xaz,"Sigmoid.Low" = yaz.low,
                    "Sigmoid.Mid" = yaz.mid,"Sigmoid.High" = yaz.high,
                    "Sample" = as.character(df$Sample)[1])) ;
}
plot_sigmoidal_fits = function(sample, bl, cper, mcount,save_eps=T,fits=T, width ,height){
  if(fits == T){
    set.fit.param.b = sigmoid_fit(sample,bin.lvl = bl, cells.per = cper, max.count = mcount, look="Either",setb = F)[1,1]
    eth.fit = sigmoid_fit(sample, bin.lvl = bl, cells.per = cper, max.count = mcount, look="Either",setb = F )
    lc.fit = sigmoid_fit(sample, bin.lvl = bl, cells.per = cper, max.count = mcount, look="LC",setb = T, param.b = set.fit.param.b )
    hc.fit = sigmoid_fit(sample, bin.lvl = bl, cells.per = cper, max.count = mcount, look="HC" ,setb = T, param.b = set.fit.param.b )
    both.fit = sigmoid_fit(sample, bin.lvl = bl, cells.per = cper, max.count = mcount, look="Both",setb = T, param.b = set.fit.param.b )
    
    eth.sig = sigmoid_fit_produce(eth.fit,max.value = 5000, step.size = 100,bval = set.fit.param.b)   
    lc.sig = sigmoid_fit_produce(lc.fit,max.value = 5000, step.size = 100,bval = set.fit.param.b)   
    hc.sig = sigmoid_fit_produce(hc.fit,max.value = 5000, step.size = 100,bval = set.fit.param.b)   
    both.sig = sigmoid_fit_produce(both.fit,max.value = 5000, step.size = 100,bval = set.fit.param.b)   
    
    
    either.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "Either")
    lc.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "LC")
    hc.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "HC")
    both.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "Both")
    
    
    if(save_eps){pdf(file=paste0(path.figs.home,sample,".sigmoid.Bcell.eps"), 
                     width=width/2.54, height=height/2.54, paper="special", bg="white",
                     fonts="Helvetica",colormode="cmyk",pointsize=1)}
    print(
      ggplot()+
        geom_line(data = lc.sig,aes(BinExp,Sigmoid.Mid),col="#E69F00")+
        geom_point(data = lc.points, aes(BinExp,propCell),col="#E69F00",size=1)+
        geom_line(data = hc.sig,aes(BinExp,Sigmoid.Mid),col="#56B4E9")+
        geom_point(data = hc.points, aes(BinExp,propCell),col="#56B4E9",size=1)+
        geom_line(data = eth.sig,aes(BinExp,Sigmoid.Mid),col="darkgrey")+
        geom_point(data = either.points, aes(BinExp,propCell),col="darkgrey",size=1)+
        geom_line(data = both.sig,aes(BinExp,Sigmoid.Mid),col="red")+
        geom_point(data = both.points, aes(BinExp,propCell),col="red",size=1)+
        theme(legend.position="none")+theme_bw()+
        scale_x_continuous(limits = c(0,mcount), breaks = c(0, 1500, 3000))+
        scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1), labels = scales::percent)+
        ylab("B cells") + xlab("Unique transcripts") +
        #scale_color_manual(values=cbpal)+
        theme(axis.title.x=element_text(family="Helvetica", size = 8),
              axis.title.y=element_text(family="Helvetica", size = 8),
              axis.text=element_text(family="Helvetica",size = 8),
              plot.title=element_text(family="Helvetica",size = 8))
      
    )
    if(save_eps){dev.off()}
  }else{
    lc.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "LC")
    hc.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "HC")
    either.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "Either")
    both.points = binCells(sam = sample,binlevel = bl, cells.per = cper,look = "Both")
    if(save_eps){pdf(file=paste0("/home/psb84/",sample,".sigmoid.Bcell.eps"), 
                     width=width/2.54, height=height/2.54, paper="special", bg="white",
                     fonts="Helvetica",colormode="cmyk",pointsize=1)}
    print(
      ggplot()+
        geom_point(data = lc.points, aes(BinExp,propCell),col="#E69F00",size=1)+
        geom_point(data = hc.points, aes(BinExp,propCell),col="#56B4E9",size=1)+
        geom_point(data = either.points, aes(BinExp,propCell),col="darkgrey",size=1)+
        geom_point(data = both.points, aes(BinExp,propCell),col="red",size=1)+
        theme(legend.position="none")+theme_bw()+
        ylab("Fraction IG Bcells") + xlab("Host mRNA, binned") +
        #scale_color_manual(values=cbpal)+
        scale_x_continuous(limits = c(0,mcount), breaks = c(0, 1500, 3000))+
        scale_y_continuous(limits = c(-0.01,1.01), breaks = c(0, 0.5, 1), labels = scales::percent)+
        theme(axis.title.x=element_text(family="Helvetica", size = 8),
              axis.title.y=element_text(family="Helvetica", size = 8),
              axis.text=element_text(family="Helvetica",size = 8),
              plot.title=element_text(family="Helvetica",size = 8))
      
    )
    if(save_eps){dev.off()}
  }
}



heavylight = function(subframe, cell){
  light.chain = subframe[subframe$Cell == cell & 
                           (subframe$Name %in% grep("^IGL",subframe$Name,value = T) |
                              subframe$Name %in% grep("^IGK",subframe$Name,value = T)),]
  
  heavy.chain = subframe[subframe$Cell == cell & 
                           subframe$Name %in% grep("^IGH",subframe$Name,value = T),]
  
  high.light.chain = as.character(as.matrix(light.chain[order(as.numeric(as.matrix(light.chain$Score)),decreasing = T),][1,1]))
  high.heavy.chain = as.character(as.matrix(heavy.chain[order(as.numeric(as.matrix(heavy.chain$Score)),decreasing = T),][1,1]))
  
  return(data.frame("Heavy" = high.heavy.chain, "Light" = high.light.chain, "Cell" = cell))
}

highhit_byregion = function(dframe, region, paired_only=T){
  subframe = data.frame(cbind(as.matrix(dframe[[paste0("best",region,"Gene")]]),
                              as.matrix(dframe[[paste0("best",region,"HitScore")]]),
                              as.matrix(dframe[["Cell"]])))
  colnames(subframe) = c("Name","Score","Cell")
  
  
  cellist = unique(as.character(as.matrix(subframe$Cell)))
  new.df = data.frame(matrix(unlist(sapply(cellist, function(cl) heavylight(subframe = subframe , cell = cl))),ncol=3,byrow = T))
  
  colnames(new.df) = c("HeavyFull","LightFull","Cell")  
  
  if (paired_only ==T){
    
    final.df = new.df[!(is.na(new.df$HeavyFull) | is.na(new.df$LightFull)),]
    final.df$Heavy = matrix(unlist(strsplit(as.character(as.matrix(final.df$HeavyFull)),split = "-")),ncol = 2,byrow = T)[,1]
    final.df$Light = matrix(unlist(strsplit(as.character(as.matrix(final.df$LightFull)),split = "-")),ncol = 2,byrow = T)[,1]
  } else{
    final.df = new.df
  }
  return(final.df)  
}

call_alignments = function(sample){
  read.cell.list = data.frame(fread(paste0("results/",sample,"/",sample,".tagged.Bcell.numerated.cell.list"),sep = "\t", header = F))
  colnames(read.cell.list) = c("readId","Cell")
  IG.stats = data.frame(fread(paste0("results/",sample,"/",sample,".Bcell.IGstats.txt"),sep = "\t", header = T))
  alignments = data.frame(fread(paste0("results/",sample,"/",sample,".tagged.Bcell.alignments.txt"),sep = "\t", header = T))
  
  
  cell.either = as.character(as.matrix(IG.stats[IG.stats$play %in% c("HC","LC","Both"),]$Cell))
  read.cell.sub.list = read.cell.list[read.cell.list$Cell %in% cell.either,]
  
  all.aligns = merge(read.cell.sub.list, alignments, by = "readId", all.x = T)
  return(all.aligns)
}

get_umi = function(sample){
  igs = data.frame(fread(paste0("results/",sample,"/",sample,".Bcell.IGstats.txt"),sep = "\t", header = T))
  exp.matrix = data.frame(fread(paste0("results/",sample,"/",sample,"_expression_matrix.txt"),sep = "\t", header = T)[,-1])
  sub.mat = data.matrix(exp.matrix[,colnames(exp.matrix) %in% unique(as.character(as.matrix(igs$Cell)))])
  return(sum(as.numeric(as.matrix(colSums(sub.mat)))))
}
