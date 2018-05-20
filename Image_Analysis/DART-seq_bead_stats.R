###
#  Authors: P. Burnham, I. De Vlaminck (2018) 
#  
#  File: DART-seq_bead_stats.R
#
#  Purpose: Use image files from accompanying python script to produce 
#  statistics file (mean intensity of beads, number of beads, etc...). 
###

#Dependencies------------------

# Packages
library(ggplot2);library(scales);

# Input variables
channel = "cy5" #Fluorescent channel used on microscope

# Paths

#Functions---------------------

get_bead_stats = function(dil, remove_extremes =F ,df = df.imagestats){
  tmp.df = df[(df$Dil == dil), ]   #recover only dilution
  
  # remove top three and bottom three points
  sort.df = tmp.df[order(tmp.df$Flu),]
  if(remove_extremes ==T){
    clip.df = sort.df[c(-1,-2,-3,-1*(nrow(sort.df)-2), -1*(nrow(sort.df)-1), -1*nrow(sort.df)),]
  }else{
    clip.df = sort.df
  }
  clip.df.test = t.test(clip.df$Flu,conf.level = 0.95)  
  final.df = data.frame("uFLU" =  clip.df.test$estimate, 
                        "maxFLU" =min(255,clip.df.test$conf.int[2]),
                        "minFLU" =max(0,clip.df.test$conf.int[1]),
                         "pmoles" = 200/dil,
                        "Number" = nrow(clip.df))
  return(final.df)
}

#Analysis-----------------------

# find all bead fluorescence property files
file.list = list.files(pattern = paste0(channel,".*.bead_details.csv"))

# create a dataframe for all bead image stats
image.list = c()

for (i in file.list){
  if(file.size(i)>0){
    image.list = rbind(image.list,read.table(i, sep=",", header = F)) ; 
  }
}

df.imagestats = data.frame(image.list) ;
colnames(df.imagestats) = c("Flu","Area","Dil","Rep") ;

# collect data for all different dilutions of toeholds

total.df = data.frame(matrix(unlist(sapply(c(0,10,15,25,50,100,200,500), 
                function(x) get_bead_stats(dil= x , remove_extremes = F , df = df.imagestats))),
                             ncol = 5, byrow = T))
colnames(total.df) = c("uFLU","maxFLU","minFLU","Mol","Number") ;
total.df[ total.df$Mol == Inf,]$Mol =0 ;

#File output------------------------------

# Table of beads statistics
write.table(total.df,"bead_imaging.stats.txt", sep="\t",quote = F, row.names = F, col.names = T)

# Figure 1b in paper
pdf(file="bead_toehold.imagstats.eps", width=5.28/2.54, 
    height=4/2.54, paper="special", bg="white",fonts="Helvetica",
    colormode="cmyk",pointsize=4)

ggplot(total.df, aes(Mol, uFLU/255)) +
  geom_errorbar(aes(ymin = minFLU/255, ymax = maxFLU/255), width = 1,alpha=0.6)+
  #geom_line() +
  geom_point(size=.75) +
  ylab("Relative intensity")+
  xlab("Toehold probes (pmoles)")+
  theme_bw()+ 
  scale_y_continuous(labels = scales::percent, limits = c(-.01,1.01), breaks = c(0,.5,1))+
  scale_x_continuous(breaks = c(0,10,20), limits = c(-4,24))+
  theme(plot.background=element_blank()) +
  theme(legend.position="none")+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

dev.off()

