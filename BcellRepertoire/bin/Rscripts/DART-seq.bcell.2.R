### Burnham, De Vlaminck: 2018
# Description: Compare IG results among the different datasets.

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")

args = commandArgs(trailingOnly = TRUE) 
sample = as.character(args[1])
umi.lb = as.numeric(args[2])
save_eps= as.logical(args[3])

# Execution / Figures ---------------------------------------------

write.table(x = quick_IG_stats(sample,umi.lb),
            file = paste0(path.stats.home, sample,".IGfractions.tab"),
            quote = F, row.names = F, col.names = F, sep="\t")

# sigmoid plot for heavy and light chains
if(save_eps){pdf(file=paste0(path.figs.home,sample,".DART-seq.Fig3D.eps"), 
                 width=4.25/2.54, height=2.8/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
plot_sigmoidal_fits(sample =sample,bl = 300,cper = 10, mcount = 3000,
                    save_eps=F,width = 5,height = 3.6 )
if(save_eps){dev.off()}
