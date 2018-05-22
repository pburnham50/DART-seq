
#####
#initialize libraries

library(Seurat) ; library(dplyr) ; library(Matrix) ; library(data.table)

#####
# Settings

args = commandArgs(trailingOnly = TRUE)
path = "/workdir/psb84/Drop_BCells/" 
path.results = paste0(path, "V5/")
figs.path = "/home/psb84/"
sam = "Bcell_1-25"
save_eps = F

#####

call_vars_Bcell = function(sam,region="C",subtypes=T){
  path.sam = paste0(path.results,sam,"/") ;
  df.reg = data.frame(fread(paste0(path.sam,sam,".Bcell.",region,"genecount.txt"))) ;
  colnames(df.reg) = c("Name","Count")
  blank =as.numeric(as.matrix(df.reg[1,2])) ; 
  df.reg = df.reg[-1,] ;
  df.reg$Fraction = df.reg$Count / sum( df.reg$Count) ;
  df.reg$Blank = blank ;
  df.reg$Region = region ;
  df.reg$Sample = sam
  if(subtypes ==F){
      df.reg$TwoCo = substr(df.reg$Name,start = 3,stop =4)
      df.agg = aggregate(Count~(TwoCo+Blank+Region+Sample),df.reg,FUN=sum)
      df.agg$Fraction = df.agg$Count/sum(df.agg$Count)
      colnames(df.agg)[1] = "Name" ;
      return(df.agg)
  }else{
      return(df.reg)
  }
}

df.ts = rbind(call_vars_Bcell("Bcell_1-100","C",subtypes = F),call_vars_Bcell("BGX3_PBMC-Reg","C",subtypes = F))
df.ts$Maj = substr(df.ts$Name,1,1)
if(save_eps){pdf(file=paste0(figs.path,"constant.pie.eps"), 
                 width=7.5/2.54, height=4.25/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}


ggplot(df.ts,aes(x="",y=Fraction,fill=Name))+
  facet_wrap(~Sample,ncol = 2)+
  geom_bar(width = 1, stat = "identity",col="black")+
  coord_polar("y",start = 0)+
  theme_bw()+ #theme(legend.position="none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        plot.title=element_text(family="Helvetica",size = 8),
        legend.text = element_text(family="Helvetica",size = 8))
if(save_eps){dev.off()}



if(save_eps){pdf(file=paste0(figs.path,"constant.bar.eps"), 
                 width=6.75/2.54, height=2.4/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}


ggplot(df.ts,aes(x=Name,y=Count,fill=Maj))+
  facet_wrap(~Sample,ncol = 2)+
  geom_bar(width = 1, stat = "identity",col="black")+
  theme_bw()+ theme(legend.position="none")+
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#D55E00"))+
  theme(axis.text.x = element_text(family="Helvetica", size = 6),
        axis.title.x = element_blank(),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        #panel.border = element_blank(),
        # panel.grid = element_blank(),
        strip.text =  element_blank(),
        strip.background = element_blank())
if(save_eps){dev.off()}


df.ts = rbind(call_vars_Bcell("Bcell_1-100","V",subtypes = T),call_vars_Bcell("BGX3_PBMC-Reg","V",subtypes = T))
df.ts$minorgroup = sapply(strsplit(df.ts$Name,"-"), `[`, 1)
df.ts$Mgroup = substr(df.ts$minorgroup,1,3)

df.ts.agg = aggregate(Count~(minorgroup+Sample),df.ts,FUN=sum)
df.ts.agg$Mgroup = substr(df.ts.agg$minorgroup,1,3)

if(save_eps){pdf(file=paste0(figs.path,sam,".variable.bar.eps"), 
                 width=6.75/2.54, height=2.75/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(df.ts,aes(x=Name,y=Count,fill=Mgroup))+
  facet_grid(Sample~Mgroup, scales ="free_x",space ="free_x")+
  geom_bar(width = 1, stat = "identity",col="black")+
  theme_bw()+ theme(legend.position="none")+
 scale_fill_manual(values = c("#56B4E9","#E69F00", "#D55E00"))+
  #scale_y_continuous(limits = c(-1,110),breaks = c(0,40,80))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        #panel.border = element_blank(),
       # panel.grid = element_blank(),
        strip.text =  element_blank(),
        strip.background = element_blank())
if(save_eps){dev.off()}



df.ts = rbind(call_vars_Bcell("Bcell_1-100","D",subtypes = T),call_vars_Bcell("BGX3_PBMC-Reg","D",subtypes = T))

if(save_eps){pdf(file=paste0(figs.path,sam,".variable.bar.eps"), 
                 width=5.6/2.54, height=5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(df.ts,aes(x=Name,y=Count,fill=Name))+
  facet_wrap(~Sample,ncol = 2)+
  geom_bar(width = 1, stat = "identity",col="black")+
  theme_bw()+ theme(legend.position="none")+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        #panel.border = element_blank(),
        # panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank())
if(save_eps){dev.off()}

