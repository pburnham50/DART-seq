library(ggplot2); library(data.table)

setwd("/workdir/psb84/Drop_BCells")

########################################
# Functions


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
  read.cell.list = data.frame(fread(paste0("V5/",sample,"/",sample,".tagged.Bcell.numerated.cell.list"),sep = "\t", header = F))
  colnames(read.cell.list) = c("readId","Cell")
  IG.stats = data.frame(fread(paste0("V5/",sample,"/",sample,".Bcell.IGstats.txt"),sep = "\t", header = T))
  alignments = data.frame(fread(paste0("V5/",sample,"/",sample,".tagged.Bcell.alignments.txt"),sep = "\t", header = T))
  
  
  cell.either = as.character(as.matrix(IG.stats[IG.stats$play %in% c("HC","LC","Both"),]$Cell))
  read.cell.sub.list = read.cell.list[read.cell.list$Cell %in% cell.either,]
  
  all.aligns = merge(read.cell.sub.list, alignments, by = "readId", all.x = T)
  return(all.aligns)
}

get_umi = function(sample){
  igs = data.frame(fread(paste0("V5/",sample,"/",sample,".Bcell.IGstats.txt"),sep = "\t", header = T))
  exp.matrix = data.frame(fread(paste0("V5/",sample,"/",sample,"_expression_matrix.txt"),sep = "\t", header = T)[,-1])
  sub.mat = data.matrix(exp.matrix[,colnames(exp.matrix) %in% unique(as.character(as.matrix(igs$Cell)))])
  return(sum(as.numeric(as.matrix(colSums(sub.mat)))))
}
########################################


hhV_regular = data.frame(highhit_byregion(call_alignments("BGX3_PBMC-Reg"),"V",paired_only = F), "Sample" = "BGX3_PBMC-Reg")

hhV_DART= data.frame(highhit_byregion(call_alignments("Bcell_1-100"),"V",paired_only = F), "Sample" = "Bcell_1-100")

hhr = hhV_DART[!is.na(hhV_DART$HeavyFull) & !is.na(hhV_DART$LightFull), ]

hhr$Count = 1

agg.hhr = aggregate(Count~(HeavyFull+LightFull), hhr, FUN=sum)
agg.hhr$LightFam = substr(x = agg.hhr$LightFull,3,3)
agg.hhr$HeavyNumber = substr(x = agg.hhr$HeavyFull,5,5)
agg.hhr$LightNumber = substr(x = agg.hhr$LightFull,5,5)

agg.agg.hhr = aggregate(Count~(HeavyNumber+LightNumber+LightFam),agg.hhr,FUN=sum)


if(save_eps){pdf(file=paste0(path.figs.home,"VarMatiching.Bcell_1to100.3.eps"), 
                 width=4.1/2.54, height=2.8/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(agg.agg.hhr,aes(LightNumber,HeavyNumber)) + 
  facet_grid(.~LightFam, space = "free",scales = "free")+geom_tile(aes(fill=Count),col="black")+
  scale_fill_gradient2(low="white",mid = "lightblue",high = "purple",midpoint = 2)+
  xlab("") + ylab("IGHV Family")+
  theme(legend.position="none")+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_blank())
if(save_eps){dev.off()}


####################

DART.1 = data.frame("Var"=c(as.character(as.matrix(hhV_DART[,1])),as.character(as.matrix(hhV_DART[,2]))),"Name"="DART" , "Count"=1)
DART.agg = aggregate(Count~(Var+Name),DART.1,FUN = sum)
Drop.1 = data.frame("Var"=c(as.character(as.matrix(hhV_regular[,1])),as.character(as.matrix(hhV_regular[,2]))),"Name"="Drop" , "Count"=1)
Drop.agg = aggregate(Count~(Var+Name),Drop.1,FUN = sum)

DART.agg$UMI = get_umi("Bcell_1-100")
Drop.agg$UMI = get_umi("BGX3_PBMC-Reg")

df = rbind(DART.agg, Drop.agg)
df$NormCount = (df$Count * 10**6)/df$UMI
df$Family = substr(df$Var,1,3)


if(save_eps){pdf(file=paste0(path.figs.home,"varsubgroups.2.eps"), 
                 width=4.8/2.54, height=2.9/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(df, aes(Var,NormCount)) +
  facet_grid(Name~Family, scales = "free_x", space = "free_x") +
  geom_bar(stat = "identity",aes(fill = Family))+
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#D55E00"))+
  ylab("Frequency")+
  xlab("Variable subtypes")+
  scale_y_continuous(breaks = c(0,50,100))+theme_bw()+  theme(legend.position="none")+ 
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text.x=element_blank(),
        axis.text.y=element_text(family="Helvetica",size = 6),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border=element_rect(colour="black",size=0.25),
        panel.background=element_rect(fill="white"),
        axis.ticks.length=unit(.05, "cm"))
if(save_eps){dev.off()}
