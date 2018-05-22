### Burnham, De Vlaminck: 2018
# Description: Used to look at variable region diversity be
# HC and LC reads.

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")
args = commandArgs(trailingOnly = TRUE)
sam.prime = as.character(args[1])
sam.compare = as.character(args[2])
save_eps = as.logical(args[3])

# Execution ---------------------------------------------

# count all variable subtypes for pairing
hhV_prime= data.frame(highhit_byregion(call_alignments(sam.prime),"V",paired_only = F), "Sample" = sam.prime)
hhV_compare = data.frame(highhit_byregion(call_alignments(sam.compare),"V",paired_only = F), "Sample" = sam.compare)


prime.1 = data.frame("Var"=c(as.character(as.matrix(hhV_prime[,1])),as.character(as.matrix(hhV_prime[,2]))),"Name"="A", "Count"=1)
prime.agg = aggregate(Count~(Var+Name),prime.1,FUN = sum)
comp.1 = data.frame("Var"=c(as.character(as.matrix(hhV_compare[,1])),as.character(as.matrix(hhV_compare[,2]))),"Name"="B" , "Count"=1)
comp.agg = aggregate(Count~(Var+Name),comp.1,FUN = sum)

prime.agg$UMI = get_umi(sam.prime)
comp.agg$UMI = get_umi(sam.compare)

df = rbind(prime.agg, comp.agg)
df$NormCount = (df$Count * 10**6)/df$UMI
df$Family = substr(df$Var,1,3)


# look at pairing in the variable region
hhr = hhV_prime[!is.na(hhV_prime$HeavyFull) & !is.na(hhV_prime$LightFull), ] # only look at cases with either HC or LC

hhr$Count = 1 # initialize count for aggregation 

agg.hhr = aggregate(Count~(HeavyFull+LightFull), hhr, FUN=sum)
agg.hhr$LightFam = substr(x = agg.hhr$LightFull,3,3)
agg.hhr$HeavyNumber = substr(x = agg.hhr$HeavyFull,5,5)
agg.hhr$LightNumber = substr(x = agg.hhr$LightFull,5,5)

agg.agg.hhr = aggregate(Count~(HeavyNumber+LightNumber+LightFam),agg.hhr,FUN=sum)

# Figures -------------------------------------------------------

if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig3E.eps"), 
                 width=4.8/2.54, height=2.9/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(df, aes(Var,NormCount)) +
  facet_grid(Name~Family, scales = "free_x", space = "free_x") +
  geom_bar(stat = "identity",aes(fill = Family))+
  scale_fill_manual(values = c("#56B4E9","#E69F00", "#D55E00"))+
  ylab("Frequency")+xlab("Variable subtypes")+
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

if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig3F.eps"), 
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