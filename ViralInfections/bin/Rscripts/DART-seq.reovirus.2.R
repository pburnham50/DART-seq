### Burnham, De Vlaminck: 2018
# Description: Used to calculate highly mutagenic sights and
# analyze transition and transversion rates in single cells.

# Settings --------------------------
source("bin/Rscripts/DARTseq.dependencies_load.R")

args = commandArgs(trailingOnly = TRUE) 
sam = as.character(args[1]) #"Lcell_T3Dpos_DARTseq_1"
path.sam = paste0(path.results, sam,"/")
min.depth = as.numeric(args[2]) #default=50, minimum depth at each SNP in the analysis
ma.low = as.numeric(args[3]) # default=0.10, min. minor allele fraction
segment = as.character(args[4]) # default="T3D_S2"
save_eps=as.logical(args[5])

# Execution -------------------------
dfol = determine_majorallele(virusSNP_qc_cell(sample = sam, cell.viability = T))

#filter SNPs based on depth of sequencing and MAF
virus.passing = nrow(dfol[dfol$TotalDepth >= min.depth & dfol$Segment == segment,])
virus.snp = nrow(dfol[dfol$TotalDepth >= min.depth &  dfol$MAF <= (1-ma.low)&dfol$Segment == segment,])
quasi.candidates.mat = dfol[dfol$TotalDepth >= min.depth &  dfol$MAF <= (1-ma.low)&dfol$Segment == segment,]

# create aggregated matrix
amlec = data.frame(rbind(as.matrix(data.frame(aggregate(A~Reference,quasi.candidates.mat,FUN = sum),"A")),
              as.matrix(data.frame(aggregate(C~Reference,quasi.candidates.mat,FUN = sum),"C")),
              as.matrix(data.frame(aggregate(G~Reference,quasi.candidates.mat,FUN = sum),"G")),
              as.matrix(data.frame(aggregate(U~Reference,quasi.candidates.mat,FUN = sum),"U"))))
colnames(amlec) = c("Reference","Count","Mutant")


agg.amlec = aggregate(as.numeric(as.matrix(Count))~Reference, amlec,FUN=sum)
colnames(agg.amlec)[2] = "Total"
nlec = merge(amlec,agg.amlec, by = "Reference")
nlec$Fraction = as.numeric(as.matrix(nlec$Count)) / as.numeric(as.matrix(nlec$Total))
nlec$newRef = factor(nlec$Reference, levels = c("U","G","C","A"))

# remove redudant cases of A>A, C>C, etc...
nle = nlec[nlec$Reference != nlec$Mutant,]
nle$mutsum = as.numeric(as.matrix(nle$Count)) / sum(as.numeric(as.matrix(nle$Count)))


# observe mutation on single cell level, based on sites described above.

qs = virusSNP_qc_cell(sample = sam, cell.viability = T)
quasi.candidates = quasi.candidates.mat$Position
qs = qs[qs$Segment == segment & qs$Position %in% quasi.candidates,]

qs.agg.cell.A = aggregate(A~(Cell+Reference),qs,FUN=sum)
qs.agg.cell.C = aggregate(C~(Cell+Reference),qs,FUN=sum)
qs.agg.cell.G = aggregate(G~(Cell+Reference),qs,FUN=sum)
qs.agg.cell.U = aggregate(U~(Cell+Reference),qs,FUN=sum)

# create aggregate matrix
all.quasi = rbind(as.matrix(data.frame(qs.agg.cell.A,"A")),
                  as.matrix(data.frame(qs.agg.cell.C,"C")),
                  as.matrix(data.frame(qs.agg.cell.G,"G")),
                  as.matrix(data.frame(qs.agg.cell.U,"U")))

all.quasi[all.quasi[,2] == "T",2] = "U"
alquasi = data.frame(all.quasi)
colnames(alquasi) = c("Cell","Ref","Count","Mut")
alquasi$CREF = paste0(alquasi$Cell, alquasi$Ref)

agg.quasi = aggregate(as.numeric(as.matrix(Count))~CREF,alquasi,FUN=sum)
colnames(agg.quasi)[2] = "Total"
alq = merge(alquasi, agg.quasi, by = "CREF")[,-1]
alq$Frac = as.numeric(as.matrix(alq$Count)) / alq$Total # conversions at s.c. basis
alq.sub = alq[alq$Total >= 20,] # only include cells with at least 20 reads

# Tables ---------------------------------------------
write.table(x = quasi.candidates.mat, file = paste0(path.stats.home,sam,".quasi_candidates.tab"),
            sep = "\t", quote = F, row.names = F, col.names = T) ;
write.table(x = nle, file = paste0(path.stats.home,sam,".quasi_conversions.bulk.tab"),
            sep = "\t", quote = F, row.names = F, col.names = T) ;
write.table(x = alq.sub, file = paste0(path.stats.home,sam,".quasi_conversions.sc.tab"),
            sep = "\t", quote = F, row.names = F, col.names = T) ;

# Figures ---------------------------------------------
# plot heat map showing frequencing of M > N
if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2F_top.eps"), 
                 width=2.75/2.54, height=2.75/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=.01)}

ggplot(nle,aes(Mutant,newRef)) + 
  geom_tile(aes(fill = mutsum*100), col = "black")+
  scale_fill_gradient(low = "yellow",high = "red")+
  scale_size(range = c(0,1))+
  xlab("Mutant")+ylab("Reference")+
  theme(legend.position = "none",
        axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8))

if(save_eps){dev.off()}

# plot the distrbution of G>N frequencies at single cell level
if(save_eps){pdf(file=paste0(path.figs.home,"DART-seq.Fig2F_bottom.eps"),
                 width=3.5/2.54, height=2.54/2.54, paper="special", bg="white",
                 fonts="Helvetica",colormode="cmyk",pointsize=1)}

ggplot(alq.sub[alq.sub$Ref == "G",],aes(Frac))+
  geom_freqpoly(aes(col=Mut),binwidth=0.05)+
  theme_bw()+theme(legend.position="none")+
  xlab("Nucleotide ratio") + ylab("Cell count")+
  scale_x_continuous(breaks = c(0,0.5,1.0),limits=c(0,1.05),labels = scales::percent)+
  scale_y_continuous(breaks=c(0,100,200))+
  theme(axis.title.x=element_text(family="Helvetica", size = 8),
        axis.title.y=element_text(family="Helvetica", size = 8),
        axis.text=element_text(family="Helvetica",size = 8),
        plot.title=element_text(family="Helvetica",size = 8),
        strip.background = element_blank(),
        strip.text = element_blank())

if(save_eps){dev.off()}