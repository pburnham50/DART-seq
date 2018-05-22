#####
#initialize libraries

library(Seurat) ; library(dplyr) ; library(Matrix) ; library(data.table)

#####
# Settings

args = c("Bcell_1-100")#commandArgs(trailingOnly = TRUE)
main.path = "/workdir/psb84/Drop_BCells/" 
figs.path = "/home/psb84/"
sam=as.character(args[1]) #sam = "9203_7858_68533_HHTNJBGX5_REG_TCGCCTTA"
sam.path = paste0(main.path,"V5/",sam,"/")
save_eps = F

#####
# Inputs

tsne.data = data.frame(fread(paste0(sam.path,sam,".allcells.descriptions.txt"), sep="\t", header = T))
ig.data = data.frame(fread(paste0(sam.path,sam,".Bcell.IGstats.txt"), sep="\t", header = T))

full.picture = merge(tsne.data, ig.data, by="Cell", all.x =T)

facs = factor(full.picture$play)
levels(facs) = c(levels(facs),"Both")
full.picture$play2 = factor(full.picture$play, levels = levels(facs))

if(save_eps){pdf(file=paste0(figs.path,sam,".HCLC_color.eps"), 
                 width=5.6/2.54, height=5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot()+
  geom_point(data = subset.data.frame(x = full.picture,
                                      subset = (Bcell==0)),
             aes(x = tSNE_1, y = tSNE_2), col= "grey",shape = 18, alpha = 0.25,size=0.25)+
  geom_point(data = subset.data.frame(x = full.picture,
                                      subset = ((Bcell==1)&(play2 == "Neither"))),
             aes(x = tSNE_1, y = tSNE_2), col = "black", alpha = 0.35,size=0.4) +
  geom_point(data = subset.data.frame(x = full.picture,
                                      subset = (Bcell==1)&(play2 != "Neither")),
             aes(x = tSNE_1, y = tSNE_2, col = play), alpha = 0.5,size=0.4)+
  theme_bw()+ theme(legend.position="none")+
  xlab("Dim. 1") + ylab("Dim. 2") +
  scale_color_manual(values=c("Both" = "red",
                              "HC" = "#56B4E9",
                              "LC" = "#E69F00"))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}


######################################################################################################

ig.data = data.frame(fread(paste0(sam.path,sam,".allcell.IGstats.txt"), sep="\t", header = T))

full.picture = merge(tsne.data, ig.data, by="Cell", all.x =T)

facs = factor(full.picture$play)
levels(facs) = c(levels(facs),"Both")
full.picture$play2 = factor(full.picture$play, levels = levels(facs))

if(save_eps){pdf(file=paste0(figs.path,sam,".HCLC_color.allcells.eps"), 
                 width=5.6/2.54, height=5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot()+
  geom_point(data = subset.data.frame(x = full.picture,
                                      subset = play2 == "Neither"),
             aes(x = tSNE_1, y = tSNE_2), col = "grey", alpha = 0.35,size=0.3) +
  geom_point(data = subset.data.frame(x = full.picture,
                                      subset = play2 != "Neither"),
             aes(x = tSNE_1, y = tSNE_2, col = play), alpha = 0.5,size=0.3)+
  theme_bw()+ theme(legend.position="none")+
  xlab("tSNE-1") + ylab("tSNE-2") +
  scale_color_manual(values=c("Both" = "red",
                              "HC" = "#56B4E9",
                              "LC" = "#E69F00"))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}

