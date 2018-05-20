### Burnham, De Vlaminck: 2018
# Description: Used to identify enrichment due to toehold sequences on reovirus genome.

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")

args = commandArgs(trailingOnly = TRUE) 
sam.1 = as.character(args[1]) # default = "Lcell_T3Dpos_Dropseq_0"
sam.2 = as.character(args[2]) # default = "Lcell_T3Dpos_DARTseq_0"
sam.3 = as.character(args[3]) # default = "Lcell_T3Dpos_DARTseq_1"
segment = as.character(args[4]) # default = "T3D_S2"
save_eps=as.logical(args[5])

# Execution ------------------------------------------

# set regular samples depth
toehold.props = data.frame(read.table(paste0(path.ref,"DART-seq.Reovirus.Toeholdinfo.csv"), header = T, sep=","))


depth_stats = rbind(virus_depth_aggregate(sam.1),
                    virus_depth_aggregate(sam.2),
                    virus_depth_aggregate(sam.3))

order = data.frame("Order"= c(1,2,3),"Sample" = c(sam.1, sam.2, sam.3))
depth_stats = merge(order,depth_stats,by="Sample")
toehold.props = merge(order,toehold.props,by="Sample")


indi_traces = rbind(toehold_enrichment_trace(sample = sam.2,stat = "Depth.Norm.UMI",region_size = 400),
                    toehold_enrichment_trace(sample = sam.1,stat = "Depth.Norm.UMI",region_size = 400))
indi_traces.agg = aggregate(Stat~(Order+Sample),indi_traces,FUN = mean)

# Metrics ---------------------------------------------------------

exp.1 = "200 bp upstream enrichment"
sum_traces.agg = aggregate(Stat~(Sample),indi_traces[indi_traces$Order <=200,],FUN = mean)
met.1 =(sum_traces.agg[sum_traces.agg$Sample == sam.2,]$Stat)/(sum_traces.agg[sum_traces.agg$Sample == sam.1,]$Stat)

exp.2 = "400 bp upstream enrichment"
sum_traces.agg = aggregate(Stat~(Sample),indi_traces[indi_traces$Order <=400,],FUN = mean)
met.2 = (sum_traces.agg[sum_traces.agg$Sample == sam.2,]$Stat)/(sum_traces.agg[sum_traces.agg$Sample == sam.1,]$Stat)

exp.3 = "S2 total enrichment"
met.3 = sum(depth_stats[depth_stats$Segment == segment & depth_stats$Sample == sam.3,]$Depth.Norm.UMI)/
          sum(depth_stats[depth_stats$Segment == segment & depth_stats$Sample == sam.1,]$Depth.Norm.UMI)

met.tab = cbind(c(exp.1,exp.2,exp.3),c(round(met.1,3),round(met.2,3),round(met.3,3)))
write.table(x = met.tab, file = paste0(path.stats.home, "DART-seq.virus.stats.tab"),
            sep = "\t",quote = F,row.names = F,col.names = F)

# Figures ---------------------------------------------------------

if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2C.eps"), width=14.5/2.54, 
                 height=3.6/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot()+
  facet_grid(Order~Segment, scales = "free_x", space = "free_x")+
  geom_area(data=depth_stats,aes(Position,(10**6 )* Depth.Norm.UMI, fill=as.factor(Order)))+
  geom_vline(data = toehold.props[toehold.props$Order == 2,],
             aes(xintercept=Start), color = "black", linetype="dashed",size=.5)+
  scale_x_continuous(breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000,4500)) +
  scale_y_continuous(breaks = c(0,40,80), limits = c(0,100)) + 
  theme_bw()+xlab("Position")+ylab("Normalized coverage")+
  theme(legend.position="none")+
  scale_fill_manual(values = c("#E69F00",cbpal[4],cbpal[3],cbpal[2]))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.text.x = element_blank(),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border=element_rect(colour="black",size=0.25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank())

if(save_eps){dev.off()}


if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2D.eps"), width=2.9/2.54, 
                 height=3.6/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot()+
  geom_line(data = indi_traces,aes(x=Order,y=(10**6)*Stat,col=Sample,fill=Segment),alpha=0.2)+
  geom_line(data = indi_traces.agg,aes(x=Order,y=(10**6)*Stat,col=Sample))+
  scale_color_manual(values = c("red","#E69F00"))+
  theme_bw()+xlab("Toehold disp. (nt)")+
  theme(legend.position="none")+
  scale_x_reverse(breaks=c(0,200,400))+
  scale_y_continuous(breaks = c(0,20,40))+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))
if(save_eps){dev.off()}


if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2E_top.eps"), width=5/2.54, 
                 height=1.6/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot()+
  facet_grid(Sample~Segment, scales = "free_x", space = "free_x")+
  geom_area(data=depth_stats[depth_stats$Segment == segment & depth_stats$Sample == sam.1,],
            aes(Position,(10**6)*Depth.Norm.UMI, fill=Sample))+
  scale_x_continuous(limits = c(0,1350),breaks = c(0,400,800,1200)) +
  scale_y_continuous(limits = c(0,(10**6)*max(depth_stats[depth_stats$Segment == segment & 
                                                            depth_stats$Sample == sam.3,]$Depth.Norm.UMI)),
                     breaks = c(0,40,80)) + 
  theme_bw()+xlab("Position (nt)")+ylab("")+
  theme(legend.position="none")+
  scale_fill_manual(values = "#E69F00")+
  theme(axis.title.x= element_blank(),
        axis.text.x= element_blank(),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border=element_rect(colour="black",size=0.25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank())

if(save_eps){dev.off()}


if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2E_bottom.eps"), width=5/2.54, 
                 height=3.5/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}
ggplot()+
  facet_grid(Sample~Segment, scales = "free_x", space = "free_x")+
  geom_area(data=depth_stats[depth_stats$Segment == segment & depth_stats$Sample == sam.3,],
            aes(Position,(10**6 )* Depth.Norm.UMI, fill=Sample))+
  geom_vline(data = toehold.props[toehold.props$Sample == sam.3,],
             aes(xintercept=Start), color = "black", linetype="dashed",size=.5)+
  scale_x_continuous(limits = c(0,1350),breaks = c(0,400,800,1200)) +
  scale_y_continuous(breaks = c(0,40,80), limits = c(0,100)) + 
  theme_bw()+xlab("Position (nt)")+ylab("")+
  theme(legend.position="none")+
  scale_fill_manual(values = c(cbpal[3]))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(1, "mm"),
        panel.border=element_rect(colour="black",size=0.25),
        panel.background=element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank())

if(save_eps){dev.off()}