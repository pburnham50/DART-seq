### Burnham, De Vlaminck: 2018

Pipeline to generate data and figures related to viral infection in the following work:
 "Simultaneous multiplexed amplicon sequencing and transcriptome profiling in single cells" by Saikia, Burnham, et al.

# Data can be download from GEO/SRA Project ID # GSE113675

# Sample description/nomenclature:
single target on each of ten gene segments, infected cells [GEO label "Lcell_T3Dpos_DARTseq_0"]
no targets, infected cells [GEO label "Lcell_T3Dpos_Dropseq_0"]
seven targets on reovirus T3D S2 segment, infected cells [GEO label "Lcell_T3Dpos_DARTseq_1"]
single target on each of ten gene segments, infected cells [GEO label "Lcell_T3Dneg_DARTseq_0"]

# Steps to generate data and figures
1) Download data and place both reads for each sample in data/ .
2) Download mm10 transcriptome files and place in references/mm10/ . Index using STAR.
3) Create the environment by running the following in the DART-seq directory:
	"conda env create --name DARTseq --file environment.yaml"
4) Activate the enviroment:
	"source activate DARTseq"
5) Change the config.reo.yaml file to samples you want to run.
6) Run snakemake (see snakemake tutorial if you are unfamiliar).

# Downstream files:

bin/Rscripts/DART-seq.reovirus.1.R : generates Figures 2C-2E and associated stats in paper.
bin/Rscripts/DART-seq.reovirus.2.R : generates Figures 2F and associated stats in paper.
bin/Rscripts/DART-seq.reovirus.3.R : generates Figures 2G,2H and associated stats in paper.
bin/Rscripts/DART-seq.reovirus.4.R : generates Figures 2I and associated stats in paper.
