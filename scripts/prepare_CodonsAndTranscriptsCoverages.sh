#!/bin/bash

# Add the modules:
module add samtools/1.9;  #samtools
module add chip-seq/1.5.5 #chipextract
module add rpfTools/1.2.0 #bamToSga

# chip-seq has been released in "G. Ambrosini; R. Dreos; S. Kumar; P. Bucher. The ChIP-Seq tools and web server: a resource for analyzing ChIP-seq and other types of genomic data. BMC genomics. 2016. DOI : 10.1186/s12864-016-3288-8."
# rpfTools has been released in "Arpat, A. B. et al. Transcriptome-wide sites of collided ribosomes reveal principles of translational pausing. Genome Res 30, 985-999 (2020). https://doi.org:10.1101/gr.257741.119"


echo "Prepare chipCor/*.vs.StartCodon.dat & chipCor/*.vs.StopCodon.dat"
echo `date`


######################################################################
### Calculate single codon coverage
# 0. Make SGA files from BAM (predicted A-site)
## same for Disome samples (keeping only 56nt to 63nt reads)

[ ! -d mapping_data/sga ] && ( mkdir mapping_data/sga )

for SID in $(cut -f1 ROS_samples.coverage | sort -u )
do
    echo $SID
    for K in $(grep $SID ROS_samples.coverage | cut -f2 )
    do
	echo "  $K"
	INFILE="mapping_data/sam/${SID}.${K}.mouse_cDNA.sorted.bam"
	OUTFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"

        if [[ "$SID" == "PD"* ]]; then
		echo "    Disome"; curShift=45; 
        	minSize=56; maxSize=63
        	samtools view -h ${INFILE} | 
			awk -v ms=${minSize} -v Ms=${maxSize} '{ if ( substr($0,1,1) == "@" || (length($10) >= ms && length($10) <= Ms)) {print $0} }' | 
			samtools view -b |
                 	samtools sort -@ 25 -n 2> /dev/null | bamToSga -i /dev/stdin -s ${curShift} > $OUTFILE
	else
                echo "    Monosome"; curShift=15; 
		samtools sort -@ 25 -n $INFILE 2> /dev/null | bamToSga -i /dev/stdin -s ${curShift} > $OUTFILE
	fi
    done
done




## 1. make sga reference file with codon coordinates
sort -k1,1 Mus_musculus.GRCm38.100.startStop.APPRIS.bed | awk '$3 < 9999999999 {split($1, n, "\|"); if (n[1] != c){print; c = n[1]}}' > Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed
awk 'BEGIN{FS=OFS="\t"}{split($1, id, "\|"); for(i=$2; i<=$3; i+=3){print id[2], "cod", i+1, "+", 1}}' Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed > Mus_musculus.GRCm38.100.codons.APPRIS_uniq.sga


## 2. counting reads per codon (for Monsome samples)
[ ! -d codonCounts ] && ( mkdir codonCounts )

for SID in $(cut -f1 ROS_samples.coverage | grep "PF" | sort -u )
do
    echo $SID
    for K in $(grep $SID ROS_samples.coverage | grep "PF" | cut -f2 )
    do
	INFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.sga"
	REF="codonCounts/ROS.mouse_cDNA.codonCounts.dat"
	OUTFILE="codonCounts/tempCounts.dat"
	[ ! -f $REF ] && ( cp Mus_musculus.GRCm38.100.codons.APPRIS_uniq.sga $REF)
	sort -k1,1 -k3,3n $REF $INFILE |
	    chipscore -c 99999999 -A "cod" -B "cDNA" -b 0 -e 2 -t 0 > $OUTFILE
	mv $OUTFILE $REF
	echo -e "$SID\t$K" >> ROS_samples.order
    done
done

## 3. remove empty codons:
echo "*** 3. remove empty codons:"
awk '{c = 0; for(i=6; i<=NF; i++) {c+=$i}; if(c > 12){print}}' codonCounts/ROS.mouse_cDNA.codonCounts.dat > codonCounts/ROS.mouse_cDNA.codonCountsShort.dat


## 4. change counts file for downstram analyses
awk 'BEGIN{FS=OFS="\t"; while( (getline < "Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed") > 0){split($1, id, "\|"); n[id[2]] = id[1]}}{print n[$1], $1, $3, $6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' codonCounts/ROS.mouse_cDNA.codonCountsShort.dat > codonCounts/ROS.mouse_cDNA.codonCountsMod.dat



######################################################################
### countig reads per transcript
echo "Prepare readCount/ outfiles..."
[ ! -d readCounts ] && ( mkdir readCounts )

for SID in $(cut -f1 ROS_samples.coverage | sort -u )
do
    echo $SID
    for K in $(grep $SID ROS_samples.coverage | cut -f2 )
    do
        echo "  $K"
        INFILE="mapping_data/sam/${SID}.${K}.mouse_cDNA.sorted.bam" 
        OUTFILE="readCounts/${SID}.${K}.readCounts.dat"
        LOGFILE="readCounts/${SID}.${K}.readCounts.log"

        curShift=15
        if [[ "$SID" == "PD"* ]]; then  echo "Disome"; curShift=45; fi
        echo "    Shift used:"$curShift

        samtools sort -@ 25 -n $INFILE 2> /dev/null | bamToBed -i - | rpf-counts -a "Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed" -s ${curShift} > $OUTFILE 2> $LOGFILE
    done
done





######################################################################
### Metagene analysis of start / stop codons (using downsampled data)
### For Monosome samples
BEDFILE=Mus_musculus.GRCm38.100.startStop.APPRIS_uniq.bed

echo "Prepare *.vs.StartCodon.dat used to generate plots startSiteAllGenes.pdf"

## Subsample to lowest coverage
N=`grep "PF" ROS_samples.coverage | cut -f3 | sort -k1,1n | head -n1`
for SID in $(grep "PF" ROS_samples.coverage | cut -f1 | sort -u )
do
    echo $SID
    for K in $(grep $SID ROS_samples.coverage | cut -f2 )
    do
        echo "  $K"
        INFILE="mapping_data/sam/${SID}.${K}.mouse_cDNA.sorted.bam"
        OUTFILE="mapping_data/sga/${SID}.${K}.mouse_cDNA.clean.Asite.Sub.sga"
        samtools sort -@ 25 -n $INFILE 2> /dev/null |
            bamToSga -i /dev/stdin -s 15 -n $N > $OUTFILE
    done
done


## start codon:
for SID in $(cut -f1 ROS_samples.coverage | grep "PF" | sort -u )
do
    echo $SID
    for K in $(grep $SID ROS_samples.coverage | cut -f2 )
    do
        TARGET=mapping_data/sga/$SID.$K.mouse_cDNA.clean.Asite.Sub.sga
        OUT=chipCor/$SID.$K.vs.StartCodon.dat
        echo "$SID $K"
        awk 'BEGIN{FS=OFS="\t"}{split($1, id, "\|"); print id[2], "start", $2+1, "+", 1}' $BEDFILE |
            sort -k1,1 -k3,3n - $TARGET |
            chipextract -A "start" -B "cDNA" -b -100 -e 1000 -w 1 -c 9999999 > $OUT
    done                                                                                                                                                                                        
done

## stop codon:
for SID in $(cut -f1 ROS_samples.coverage | grep "PF" | sort -u )
do
    echo $SID
    for K in $(grep $SID ROS_samples.coverage | cut -f2 )
    do
        TARGET=mapping_data/sga/$SID.$K.mouse_cDNA.clean.Asite.Sub.sga
        OUT=chipCor/$SID.$K.vs.StopCodon.dat
        echo $I
        awk 'BEGIN{FS=OFS="\t"}{split($1, id, "\|"); print id[2], "stop", $3, "+", 1}' $BEDFILE |
            sort -k1,1 -k3,3n - $TARGET |
            chipextract -A "stop" -B "cDNA" -b -1000 -e 100 -w 1 -c 9999999 > $OUT
    done
done


