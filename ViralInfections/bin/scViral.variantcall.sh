#!/bin/bash
# scViral.variantcall.sh

#### Description
# used with barcoded .sam file to seperate into individual sorted bam files,
# then use bamread-count to determine variants at each position.

#### References
reoref=references/Reo/reo_T3D.fasta

#### Programs
bamreadcount=bin/programs/bam-readcount-master/bin/bam-readcount
samtools=/programs/samtools-1.6/bin/samtools

#### Inputs
sample=$1
sampath=results/$sample

#### Execute

mkdir $sampath/cells

rm -f $sampath/$sample.cell.list
cut -d' ' -f 3 $sampath/$sample.virus.barcoded.sam | sort | uniq > $sampath/$sample.cell.list

# individual cells get separated, sent to sorted bam, and read into variable files

cat $sampath/$sample.cell.list | while read line
do
	grep -w $line $sampath/$sample.virus.barcoded.sam | cut -d' ' -f 1,5- |\
		 sed 's/ /\t/g' > $sampath/cells/$sample.virus.barcoded.$line.tmp.sam ;
	cat $sampath/$sample.virusaln.samhead $sampath/cells/$sample.virus.barcoded.$line.tmp.sam | \
		samtools view -bS -q 55 -f 0 - > $sampath/cells/$sample.$line.bam ;
	$samtools sort -o $sampath/cells/$sample.$line.sort.bam $sampath/cells/$sample.$line.bam ;
	$samtools index $sampath/cells/$sample.$line.sort.bam ;
	$samtools depth -Q 55 -q 13 $sampath/cells/$sample.$line.sort.bam | \
		awk -v a=$sample -v b=$line '{print a"\t"b"\t"$0}' > $sampath/cells/$sample.$line.depth;
	$bamreadcount -q 55 -b 13 -w 0 -f $reoref $sampath/cells/$sample.$line.sort.bam | \
		sed 's/:/\t/g' | cut -f 1,2,3,20,34,48,62 | \
		awk -v a=$line '{print a"\t"$0}' > $sampath/cells/$sample.$line.vars ;
done

cat $sampath/cells/$sample.*.vars > $sampath/$sample.virvar.tab ;
cat $sampath/cells/$sample.*.depth > $sampath/$sample.virdepth.tab ;
rm -fr $sampath/cells
