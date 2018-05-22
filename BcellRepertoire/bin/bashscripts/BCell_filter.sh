#!/bin/bash

################################################################################################
# Software

MIXCR=bin/mixcr-2.1.5/mixcr

###########################################################################################################################
# enter inputs here
SAMPLE=$1
SAMDIR=V5/$SAMPLE #### directory with the data here
DATADIR=/workdir/Data/SingleCell
CELL_IDS=$SAMDIR/$SAMPLE.Bcell.barcodes.txt  # file with list of B cell IDs
INDIVIDUAL_CELLS=$SAMDIR/Bcellcodes # where to store individual cells
#CUST_PRIMERS=/workdir/fw262/BashScripts/V_C_primers.txt # file with primers
FASTQGZ_R1=$DATADIR/$SAMPLE\_R1.fastq.gz
FASTQGZ_R2=$DATADIR/$SAMPLE\_R2.fastq.gz
###########################################################################################################################
IND_FILES=$SAMDIR/$INDIVIDUAL_CELLS/*

#create the unaligned bam files for r1 and r2 separately
#java -Djava.io.tmpdir=temp -jar /programs/picard-tools-2.1.1/picard.jar FastqToSam\
#		F1=$FASTQGZ_R1 \
#		OUTPUT=$SAMDIR/unaligned_r1.bam \
#		SM=DS 
#java -Djava.io.tmpdir=temp -jar /programs/picard-tools-2.1.1/picard.jar FastqToSam\
#		F1=$FASTQGZ_R2 \
#		OUTPUT=$SAMDIR/unaligned_r2.bam \
#		SM=DS 
#
#echo "Rearranging BAM file so that the sequences are at the front of the line for read 1."

# rearrange bam file so sequences are at the front
#samtools view $SAMDIR/unaligned_r1.bam | cut -f 1-9,11,12 > $SAMDIR/tmp2.sam
#samtools view $SAMDIR/unaligned_r1.bam | cut -f 10 > $SAMDIR/tmp1.sam
#paste -d"\t" $SAMDIR/tmp1.sam $SAMDIR/tmp2.sam | sort > $SAMDIR/rearrangedSamSorted_r1.sam
#rm -f $SAMDIR/tmp1.sam
#rm -f $SAMDIR/tmp2.sam
#
#echo "Sorting non rearranged BAM file for read 2."
#
#samtools view unaligned_r2.bam | sort > sortedUnaligned_r2.sam
#
#echo "Now looking for cell IDs of interest"
#
-mkdir $INDIVIDUAL_CELLS
while read p; do
    look $p rearrangedSamSorted_r1.sam >> $INDIVIDUAL_CELLS/$p.sam
done <$CELL_IDS

# read1.sam contains read 1's for cell of interest (specified in CELL_IDS)
for f in $IND_FILES
do
	rm -f ${f}_readIDsofInterest.sam
	cut -f 2 $f > ${f}_tmp1.sam # contains readIDs
	cut -f 1 $f > ${f}_tmp2.sam # contains sequences
	paste -d"\t" ${f}_tmp1.sam ${f}_tmp2.sam > ${f}_readIDsofInterest.sam
	rm -f ${f}_tmp1.sam 
	rm -f ${f}_tmp2.sam
	rm -f $f
done # ${f}_readIDsofInterest.sam or each file in $IND_FILES should have read id in column1 and 20 base sequence as column2


echo "Looking for read IDs of interest from read2"

for f in $IND_FILES
do
	rm -f ${f}_read2ofInterest.sam
	while read x; do # look through each line of _readIDsofInterest.sam and find match in sortedUnaligned_r2.sam  		
		arr=($x) # x is a string right now		
		readIDsOnly=$(echo ${arr[0]})
		baseSeq20=$(echo ${arr[1]})	
		look $readIDsOnly sortedUnaligned_r2.sam | cut -f 1,10 > ${f}_read2ofInterest.sam
		sed "s/$/ \t ${baseSeq20}/" ${f}_read2ofInterest.sam >> ${f}_read2ofInterest_with_cellIDUMI.sam
	done < $f # f is the file
	rm -f ${f}_read2ofInterest.sam
	rm -f $f # this removes _readIDsofInterest.sam which also contains 20 base sequences as column2
done

echo "Now we have reads of read 2 separated by each cell."

echo "Creating file with reads of interest only as well as a read ID formated for fasta"
for f in $IND_FILES
do
	cut -f 2 $f > ${f}_readSeqFinal.txt
	cut -f 1 $f > ${f}_readIDsFinal.txt
	cut -f 3 $f > ${f}_cellbar_UMI.txt
	sed -i -e 's/^/>/' ${f}_readIDsFinal.txt # adds the > sign in front of the ID for fasta file
	paste -d : ${f}_readIDsFinal.txt ${f}_cellbar_UMI.txt > ${f}.tmp
	paste -d"\n" ${f}.tmp ${f}_readSeqFinal.txt > ${f}.fasta
	rm -f ${f}_readSeqFinal.txt ${f}_readIDsFinal.txt ${f}_cellbar_UMI.txt ${f}.tmp
	rm -f ${f}
done

echo "fasta files created for each cell of interest"
for f in $IND_FILES
do
    [ -f "$f" ] || continue
    mv "$f" "${f//.sam_readIDsofInterest.sam_read2ofInterest_with_cellIDUMI.sam/}" # rename the files here
done


echo "creating file with all cell reads in fasta format"
rm -f $INDIVIDUAL_CELLS/allCells.fasta
cat *.fasta > allCells.fasta_all

echo "Running mixCR to align sequences with IG genes for each individual B cell."

# run mixCR on individual cells
for f in $DATA_DIRECTORY/$INDIVIDUAL_CELLS/*.fasta
do

	$MIXCR align -f -s hsa -p rna-seq -OallowPartialAlignments=true ${f} ${f}.vdjca > ${f}.reportTxt
	
	$MIXCR assemblePartial -f ${f}.vdjca ${f}.rescued1.vdjca
	$MIXCR assemblePartial -f ${f}.rescued1.vdjca ${f}.rescued2.vdjca
	$MIXCR assemble -f ${f}.rescued2.vdjca ${f}.clns

	$MIXCR exportClones -f ${f}.clns ${f}.clones.txt
	$MIXCR exportClones -f -c TRB ${f}.clns ${f}.clones.TRB.txt
	$MIXCR exportClones -f -c IG ${f}.clns ${f}.clones.IGH.txt
	$MIXCR exportAlignments -f -c IG ${f}.rescued2.vdjca ${f}.alignmentsAll.txt
	$MIXCR exportAlignments -f -c IG -vGene -vHitScore -dGene -dHitScore -jGene -jHitScore -cGene -cHitScore ${f}.rescued2.vdjca ${f}.alignments.txt

	$MIXCR exportAlignmentsPretty ${f}.rescued2.vdjca ${f}.aligned_view.txt
done


###########################################
#################### igBlast analysis below 
######################## comment out if not necessary (probably not necessary since mixCR above should work)
#echo "Running igBLAST to align sequence with IG genes for each cell."
#$DIRREADS=$DATA_DIRECTORY/$INDIVIDUAL_CELLS

#cd /workdir/fw262/tools/ncbi-igblast-1.7.0/bin # go to folder with igblast and required databases
#IND_FILES=$DATA_DIRECTORY/$INDIVIDUAL_CELLS/*
#for f in $IND_FILES # run igblastn for each individual cell fasta file
#do
#	./igblastn -germline_db_V human_IG_V -germline_db_D human_IG_D -germline_db_J human_IG_J -organism human -query $DIRREADS/$f -auxiliary_data optional_file/human_gl.aux -out $DIRREADS/${f}_igblast.txt -outfmt '7' # run igblast here
#done

#echo "Go through each igblast result and extract information about total queries, identifiable CD3s, and unique clonal types"
#IND_FILES_IGBLAST=$DATA_DIRECTORY/$INDIVIDUAL_CELLS/*fasta_igblast.txt # find ONLY individual cell results
#rm -f igblastResultData.txt
#echo -e "Total queries" '\t' "Total identifiable CDR3" '\t' "Total unique clonotypes" >> igblastResultData.txt
#for ifg in $IND_FILES_IGBLAST # look for only files that end in _igblast.txt
#do
#if [["$ifg" -ef "${DATA_DIRECTORY}/$INDIVIDUAL_CELLS/allCells.fasta_all_igblast.txt"]] ; then
#	read -p "Press [Enter] key to start backup..."
#	all1=$(grep "Total queries =" $ifg | cut -d\  -f4) 
#	all2=$(grep "Total identifiable CDR3 =" $ifg | cut -d\  -f5) 
# 	all3=$(grep "Total unique clonotypes =" $ifg | cut -d\  -f5) 
#else
#	TMP1=$(grep "Total queries =" $ifg | cut -d\  -f4) 
#	TMP2=$(grep "Total identifiable CDR3 =" $ifg | cut -d\  -f5) 
#	TMP3=$(grep "Total unique clonotypes =" $ifg | cut -d\  -f5) 
#	echo -e $TMP1'\t'$TMP2'\t'$TMP3 >> igblastResultData.txt
#fi
#done
#echo -e $all1'\t'$all2'\t'$all3 >> igblastResultData.txt # add all Cell results at the very end


