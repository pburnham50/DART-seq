### Bcell_IG_comparison.R

#####
# Description
#
# Compare IG results among the different datasets.
# First need to move all of the "*.Bcell.IGstats.txt"
# into the folder group_comparison/ off the main path.
#####
# Libraries
save_eps=F
library(ggplot2); library(data.table); library(caTools);

impanno = function(sam){
  df = data.frame(fread(paste0("group_comparison/both/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
  df$Sample = sam ; 
  df$CellOrder = nrow(df)-rank(as.numeric(as.matrix(df$IGdensity)),ties="first",na.last =F)
  return(df)
}

tablify = function(sam){
  df = data.frame(fread(paste0("group_comparison/both/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
  tab = as.matrix(table(factor(df$play,levels = c("Both","HC","LC","Neither")))) ;
  newdf = data.frame("Inc" = rownames(tab), "Count" = as.numeric(tab)) ;
  newdf$Prop = newdf$Count / sum(newdf$Count) ; 
  newdf$Sample = sam ; 
  return(newdf)
}

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

binCells = function(sam, binlevel = 250, look = "LC"){
  df = data.frame(fread(paste0("group_comparison/both/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
  df$Binexp = cut(df$mRNA,seq(0,4000,binlevel),labels = F)*binlevel ;
  levs = seq(0,4000,binlevel) ;
  finale = sapply(levs, function(x) cell_frac(mrna.bin = x, df, look = look))
  nce = data.frame(levs, t(finale), sam, look)   
  colnames(nce) = c("BinExp","nCell","propCell", "Sample", "IGtype")
  return(nce)
}

#####
#
sample.list = as.character(matrix(unlist(strsplit(list.files("group_comparison/both/"),"\\.")),
                                  ncol=4,byrow = T)[,1])

allsams.df = impanno(sample.list[1])

for (i in sample.list[-1]){
  allsams.df = rbind(allsams.df,impanno(i)) ;
}


tab.df = tablify(sample.list[1])

for (i in sample.list[-1]){
  tab.df = rbind(tab.df,tablify(i)) ;
}




path.figs = "/home/psb84/"
save_eps=F



if(save_eps){pdf(file=paste0(path.figs,"IG.histo.eps"), width=7/2.54, height=7.5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}


ggplot(allsams.df,aes(play))+geom_histogram(stat="count",aes(fill=Sample), position = "dodge")+
  theme(plot.background=element_blank()) + theme(legend.position="none")+
   xlab("IG presence") + ylab("Number of B cells") +
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
         axis.title.y=element_text(family="Helvetica", size = 8),
         axis.text=element_text(family="Helvetica",size = 8),
         plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}


if(save_eps){pdf(file=paste0(path.figs,"IG.histo.eps"), width=7/2.54, height=7.5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}


ggplot(allsams.df,aes(Sample))+
  geom_histogram(stat="count",aes(fill=play), position = "stack")+
  theme(plot.background=element_blank()) + #theme(legend.position="none")+
  xlab("IG presence") + ylab("Number of B cells") +
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}

tab.df$Sample2 = factor(tab.df$Sample, levels = c("Bcell_1-25", "Bcell_1-100","BGX3_PBMC-Reg",
                                                  "9203_7858_68534_HHTNJBGX5_V1-1-100_CTAGTACG",
                                                "9203_7858_68535_HHTNJBGX5_V1-1-500_TTCTGCCT",
                                                "9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA"))
tab.df = tab.df[!is.na(tab.df$Sample2),]
library(scales)


if(save_eps){pdf(file=paste0(path.figs,"IG.bar.2.eps"), width=6.8/2.54, height=2.85/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot(tab.df,aes(x=Sample2,y = Prop))+
  geom_bar(stat="identity",aes(fill=Inc))+
  scale_x_discrete(labels = c("Bcell_1-25" = "DART-seq B, 8.0 pm",
                              "Bcell_1-100" = "DART-seq B, 2.0 pm",
                              "BGX3_PBMC-Reg" = "Drop-seq B",
                              "9203_7858_68534_HHTNJBGX5_V1-1-100_CTAGTACG" = "DART-seq A, 2.0 pm",
                              "9203_7858_68535_HHTNJBGX5_V1-1-500_TTCTGCCT" = "DART-seq A, 0.4 pm",
                              "9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA" = "Drop-seq A"))+
  scale_fill_manual(values = c("red","#56B4E9","#E69F00","black")) +
  xlab("") + ylab("B cell, IG capture") + scale_y_continuous(labels = scales::percent, breaks = c(0,.5,1))+
  theme(plot.background=element_blank()) +theme_bw()+theme(legend.position="none")+ coord_flip()+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(family="Helvetica", size = 8),
        axis.text.x=element_text(family="Helvetica",size = 8),
        axis.text.y=element_text(family="Helvetica",size = 6),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}

if(save_eps){pdf(file=paste0(path.figs,"IG.dodge.eps"), width=14/2.54, height=7.5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}


ggplot(tab.df,aes(x=Sample2,y = Count))+
  geom_bar(stat="identity",aes(fill=Inc),col="black",position="dodge")+
  scale_x_discrete(labels = c("9203_7858_68534_HHTNJBGX5_V1-1-100_CTAGTACG" = "DART-seq A\n2.0 pm",
                              "9203_7858_68535_HHTNJBGX5_V1-1-500_TTCTGCCT" = "DART-seq A\n0.4 pm",
                              "9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA" = "Drop-seq A",
                              "Bcell_1-25" = "DART-seq B\n8.0 pm",
                              "Bcell_1-100" = "DART-seq B\n2.0 pm",
                              "BGX3_PBMC-Reg" = "Drop-seq B"))+
  theme(plot.background=element_blank()) +#theme(legend.position="none")+
  xlab("") + ylab("Number of Bcells") + 
  scale_fill_manual(values = c("red","#56B4E9","#E69F00","black")) +
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text.x = element_blank(),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}




# 
# 
# if(save_eps){pdf(file=paste0(path.figs,"IG.histo.eps"), width=7/2.54, height=7.5/2.54, paper="special", bg="white",
#                  fonts="Helvetica",colormode="cmyk",pointsize=1)}
# 

bin.df = binCells(sample.list[1],binlevel = 200,look = "LC")

for (i in sample.list[-1]){
  bin.df = rbind(bin.df,binCells(i,binlevel = 200,look = "LC")) ;
}

cbpal =  c("black","#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7","red","darkgreen")


####################################################################################################################################


binCells = function(sam, binlevel = 250, max.count = 4000, look = "LC", cells.per = 10){
  df = data.frame(fread(paste0("group_comparison/both/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
  df$Binexp = cut(df$mRNA,seq(0,max.count,binlevel),labels = F)*binlevel ;
  levs = seq(0,max.count,binlevel) ;
  finale = sapply(levs, function(x) cell_frac(mrna.bin = x, df, look = look))
  nce = data.frame(levs, t(finale), sam, look)   
  colnames(nce) = c("BinExp","nCell","propCell", "Sample", "IGtype")
  return(nce[nce$nCell >= cells.per,])
}

sam = sample
df = data.frame(fread(paste0("group_comparison/both/",sam,".Bcell.IGstats.txt"),sep="\t",header = T))
df.sub = df[df$mRNA >= 1200, ]
df.hits = data.frame(table(df.sub$play))
sum(df.hits[df.hits$Var1 != "Neither",]$Freq)/sum(df.hits$Freq)
sum(df.hits[df.hits$Var1 == "Both",]$Freq)/sum(df.hits$Freq)

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
    
sample=sample.list[5]
  
bl = 300
cper = 8
mcount = 5000

######

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
      
    
    if(save_eps){pdf(file=paste0("/home/psb84/",sample,".sigmoid.Bcell.eps"), 
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

sample="Bcell_1-100"
if(save_eps){pdf(file=paste0("/home/psb84/",sample,".sigmoid.Bcell.eps"), 
                 width=4.25/2.54, height=2.8/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
plot_sigmoidal_fits(sample =sample,bl = 300,cper = 10, mcount = 3000,save_eps=F,width = 5,height = 3.6 )
if(save_eps){dev.off()}

sample="BGX3_PBMC-Reg"
if(save_eps){pdf(file=paste0("/home/psb84/",sample,".sigmoid.Bcell.eps"), 
                 width=4.25/2.54, height=2.8/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
plot_sigmoidal_fits(sample =sample,bl = 300,cper = 10, mcount = 3000,save_eps=F,width = 5,height = 3.6)
if(save_eps){dev.off()}

