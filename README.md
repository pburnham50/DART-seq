# DART-seq

Pipeline to generate data and figures related to viral infection in the following work:
 "Simultaneous multiplexed amplicon sequencing and transcriptome profiling in single cells" by Saikia, Burnham, et al.,
which is currently under review.

Files related to DART-seq technique including:
A) Image_Analysis - 	scripts for fluorescence quantification of modified Dropseq beads.
B) ViralInfections - 	Snakemake pipeline for host and virus transcription. Includes rules
				to identify viral fragments from T3D orthoreovirus and identify 
				mutations. Will also generate figures related to paper.
C) BcellRepertoire - 	Snakemake pipeline for host transcription in PBMCs, identification of
				B lymphocytes in PBMCs and subsequent determination of 
				heavy and light chain elements. Will also generate figures
				related to the paper.
