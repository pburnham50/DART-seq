library(dplyr); library(ggplot2)


wd = "/workdir/psb84/Drop_BCells/V5/Bcell_1-100/"
save_eps=F
a1 = data.frame(read.table(paste0(wd, "a.b"),sep=",", header =F))
a2 =  data.frame(read.table(paste0(wd, "c.b"),sep=",", header =F))

a12 = merge(a1,a2,by = "V1")
colnames(a12) = c("ReadNumber","Vhit","Vscore",
                  "Dhit","Dscore","Jhit","Jscore",
                  "Chit","Cscore","Cell")

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
  subframe = data.frame(cbind(as.matrix(dframe[[paste0(region,"hit")]]),
                        as.matrix(dframe[[paste0(region,"score")]]),
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

hhr = highhit_byregion(a12,"V")
hhr$Count = 1

agg.hhr = aggregate(Count~(Heavy+Light), hhr, FUN=sum)
agg.hhr$LightFam = substr(x = agg.hhr$Light,3,3)
agg.hhr$HeavyNumber = substr(x = agg.hhr$Heavy,5,5)
agg.hhr$LightNumber = substr(x = agg.hhr$Light,5,5)



if(save_eps){pdf(file=paste0(path.figs.home,"VarMatiching.Bcell_1to100.2.eps"), 
                 width=4.75/2.54, height=3.25/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(agg.hhr,aes(LightNumber,HeavyNumber)) + 
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


hh.V = highhit_byregion(a12,"V",F)
hh.D = highhit_byregion(a12,"D",F)
hh.J = highhit_byregion(a12,"J",F)
hh.C = highhit_byregion(a12,"C",F)

vdj = Reduce(function(x, y) merge(x, y, all=TRUE, by = "Cell"), list(hh.V, hh.D, hh.J, hh.C))

vdj.heavy = vdj[,c(1,2,4,6,8)]
colnames(vdj.heavy) = c("Cell","V","D","J","C")
vj.light = vdj[,c(1,3,7,9)]
colnames(vj.light) = c("Cell","V","J","C")


VDJ.df =data.frame("Stat" = c(nrow(vdj.heavy[!is.na(vdj.heavy$V),]),
                              nrow(vdj.heavy[!is.na(vdj.heavy$V)&!is.na(vdj.heavy$D),]),
                              nrow(vdj.heavy[!is.na(vdj.heavy$V)&!is.na(vdj.heavy$D)&!is.na(vdj.heavy$J),]),
                            #  nrow(vdj.heavy[!is.na(vdj.heavy$V)&!is.na(vdj.heavy$D)&!is.na(vdj.heavy$J)&!is.na(vdj.heavy$C),]),
                              nrow(vj.light[!is.na(vj.light$V),]),
                              nrow(vj.light[!is.na(vj.light$V)&!is.na(vj.light$J),])),
                        #      nrow(vj.light[!is.na(vj.light$V)&!is.na(vj.light$J)&!is.na(vj.light$C),])),
                   "Metric" = c("V","V+D","V+D+J","V","V+J"),
                   "Chain" = c(rep("Heavy",3),rep("Light",2)))




if(save_eps){pdf(file=paste0(path.figs.home,"VDJ.continuity.eps"), 
                 width=5/2.54, height=3/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(VDJ.df,aes(Metric,Stat))+geom_bar(stat="identity", aes(fill = Chain),col="black")+
  facet_grid(Chain~., scales = "free", space = "free")+ coord_flip()+
  scale_fill_manual(values = c("#56B4E9","#E69F00"))+
  theme(legend.position="none")+ xlab("")+ylab("Count")+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_blank())

if(save_eps){dev.off()}





###########################################


var_hits = as.matrix(table(c(as.character(as.matrix(hh.V$HeavyFull)), 
  as.character(as.matrix(hh.V$LightFull)))))

var.df = data.frame("Name"=rownames(var_hits), "Hits" = var_hits)
var.df$MG = substr(var.df$Name, start = 1, stop=3)

var.df$Name2 = factor(var.df$Name, levels = var.df[order(var.df$Hits,decreasing = T),]$Name)




if(save_eps){pdf(file=paste0(path.figs.home,"varsubgroups.eps"), 
                 width=4.4/2.54, height=2.75/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(var.df,aes(Name2,Hits/max(Hits)))+geom_bar(stat="identity", aes(fill = MG),col="black",size=.25)+
  facet_grid(~MG, scales = "free", space = "free")+ 
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#D55E00"))+
  theme(legend.position="none")+ xlab("")+ylab("Count")+
  xlab("Variable subtypes")+ylab("Normalized frequency")+
  scale_y_continuous(breaks=c(0,.5,1), limits = c(-.01,1.2))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_blank())

if(save_eps){dev.off()}


