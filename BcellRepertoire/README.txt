### Burnham, De Vlaminck: 2018

Pipeline to generate data and figures related to B cell repertoire in the following work:
 "Simultaneous multiplexed amplicon sequencing and transcriptome profiling in single cells" by Saikia, Burnham, et al.

# Data can be download from GEO/SRA Project ID # GSE113675

# Sample description/nomenclature:
Human PBMCs, standard Dropseq [GEO label "PBMC_Dropseq_0"]
Human PBMCs, targeting 6 heavy chain and 2 light chain constant regions [GEO label "PBMC_DARTseq_0"]

# Steps to generate data and figures
1) Download data and place both reads for each sample in data/ .
2) Download hg19 transcriptome files and place in references/hg19/ . Index using STAR.
3) Create the environment by running the following in the DART-seq directory:
	"conda env create --name DARTseq --file environment.yaml"
4) Activate the enviroment:
	"source activate DARTseq"
5) Change the config.bcell.yaml file to samples you want to run.
6) Run snakemake (see snakemake tutorial if you are unfamiliar).

# Downstream files:

bin/Rscripts/DART-seq.bcell.1.R : generates Figures 3C and associated stats in paper.
bin/Rscripts/DART-seq.bcell.2.R : generates Figures 3D and associated stats in paper.
bin/Rscripts/DART-seq.bcell.3.R : generates Figures 3E,3F and associated stats in paper.

# Commands to run for downstream:
"Rscript bin/Rscripts/DART-seq.bcell.1.R PBMC_DARTseq T ; "
"Rscript bin/Rscripts/DART-seq.bcell.2.R PBMC_DARTseq_0 1200 T ; "
"Rscript bin/Rscripts/DART-seq.bcell.2.R PBMC_Dropseq_0 1200 T ; "
"Rscript bin/Rscripts/DART-seq.bcell.3.R PBMC_DARTseq_0 PBMC_Dropseq_0 T ; "
