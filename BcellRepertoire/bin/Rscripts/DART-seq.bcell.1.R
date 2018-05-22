### Burnham, De Vlaminck: 2018
# Description: Used to plot a tSNE of the sample of interest
# and highlight cells carrying LC reads, HC reads, and both
# HC and LC reads.

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")

args = commandArgs(trailingOnly = TRUE)
sam=as.character(args[1]) 
save_eps = as.logical(args[2])
sam.path = paste0(path.results,sam,"/")

# Execution ---------------------------------------------
tsne.data = data.frame(fread(paste0(sam.path,sam,".allcells.descriptions.txt"), sep="\t", header = T))
ig.data = data.frame(fread(paste0(sam.path,sam,".Bcell.IGstats.txt"), sep="\t", header = T))

full.picture = merge(tsne.data, ig.data, by="Cell", all.x =T)

facs = factor(full.picture$play)
levels(facs) = c(levels(facs),"Both") # accounts for case of 0 cells with both heavy and light chains i
full.picture$play2 = factor(full.picture$play, levels = levels(facs))

# Figures -------------------------------------------------------
if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig3C.eps"), 
                 width=5.6/2.54, height=5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot()+
  geom_point(data = subset.data.frame(x = full.picture, subset = (Bcell==0)),
             aes(x = tSNE_1, y = tSNE_2), col= "grey",shape = 18, alpha = 0.25,size=0.25)+
  geom_point(data = subset.data.frame(x = full.picture, subset = ((Bcell==1)&(play2 == "Neither"))),
             aes(x = tSNE_1, y = tSNE_2), col = "black", alpha = 0.35,size=0.4) +
  geom_point(data = subset.data.frame(x = full.picture, subset = (Bcell==1)&(play2 != "Neither")),
             aes(x = tSNE_1, y = tSNE_2, col = play), alpha = 0.5,size=0.4)+
  theme_bw()+ theme(legend.position="none")+ xlab("tSNE-1") + ylab("tSNE-2") +
  scale_color_manual(values=c("Both" = cbpal[1], "HC" = cbpal[2],"LC" = cbpal[3]))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}