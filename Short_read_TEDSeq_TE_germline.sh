#!/bin/bash
####################################################################################
########################TED-SEQ Short-reads germline pipeline#######################
####################################################################################



####################################################################################
################################Required Tools######################################
####################################################################################
#bbmap (Version 38.18)
#cutadapt (Version 2.6)
#r-dplyr (Version 1.1.3)
#bedtools (Version 2.27.1)
#samtools (Version 1.20)
#bowtie2 (Version 2.20.1)

####################################################################################
###################################Settings#########################################
####################################################################################
show_help() {
    echo "TED-SEQ Short-reads germline pipeline"
    echo "Usage: bash Short_read_TEDSeq_TE_germline.sh [options]"
    echo
    echo "Options:"
    echo "  --sample        Name of the sample (required)."
    echo "  --read_1        Full path to first FASTQ read file (required)."
    echo "  --read_2        Full path to second FASTQ read file (required)."
    echo "  --refDir        Directory containing the reference genome (required)."
    echo "  --outDir        Output directory (required)."
    echo "  --te_family     Desired family analyzed (default: ATCOPIA93)."
    echo "  --cores         Number of cores to use (default: 20)."
    echo "  --min_cov       Minimum coverage for TE detection (30 for ATCOPIA93, 20 for ATENSPM3, 5 for VANDAL21)."
    echo "  --min_ratio     Minimum ratio (default: 0.8)."
    echo "  --libsize       Library size (default: 100)."
    echo "  --ref_TE        Reference TE bed file (required)."
    echo "  --scriptDir     Directory containing the scripts and the flanking sequences (required)."
    echo "  --help          Show this help message."
    echo
}

# Default values
file=""
read_1=""
read_2=""
refDir=/mnt/data2/leduque/ONT/reference/TAIR10_allChr
outDir=""
te_family="ATCOPIA93"
CORES=20
min_cov=30
min_ratio=0.8
libsize=100
ref_TE=""
scriptDir=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample) file=$2; shift 1;;
        --read_1) read_1=$2; shift 1;;
        --read_2) read_2=$2; shift 1;;
        --refDir) refDir=$2; shift 1;;
        --outDir) outDir=$2; shift 1;;
        --te_family) te_family=$2; shift 1;;
        --cores) CORES=$2; shift 1;;
        --min_cov) min_cov=$2; shift 1;;
        --min_ratio) min_ratio=$2; shift 1;;
        --libsize) libsize=$2; shift 1;;
        --ref_TE) ref_TE=$2; shift 1;;
        --scriptDir) scriptDir=$2; shift 1;;
        --help) show_help; exit 0;;
        *) echo "Unknown parameter passed: $1"; show_help; exit 1;;
    esac
    shift
done

if [[ -z "$file" || -z "$read_1" || -z "$read_2"|| -z "$refDir" || -z "$outDir"|| -z "$ref_TE"|| -z "$scriptDir" ]]; then
    echo "Error: lacking required arguments."
    exit 1
fi

if [[ -z "$outDir" ]]; then
    outDir=/mnt/data7/pvendrell/results_new_pipeline/$file
fi

mkdir -p $outDir/$file


# Display chosen options
echo "###############Running TED-SEQ Short-reads germline pipeline################"
echo ""
echo "############################Files###########################"
echo ""
echo "Sample: $file"
echo "Read 1: $read_1"
echo "Read 2: $read_2"
echo "Reference Directory: $refDir"
echo "Output Directory: $outDir/$file"
echo "TE Family: $te_family"
echo "Cores: $CORES"
echo "Minimum Coverage: $min_cov"
echo "Minimum Ratio: $min_ratio"
echo "Library Size: $libsize"
echo "Reference TE File: $ref_TE"
echo "Script Directory: $scriptDir"
echo ""
outTrim=$outDir/$file/Trimmed
mkdir  $outTrim

echo "############################Starting Pipeline############################"

####################################################################################
###############Subsampling and PCR duplicate removal and Clip TE ###################
####################################################################################


  # removal and Clip TE 
  #flanking_sequence.fa contains from the primer to almost the extremity of the desired family
  #flanking_sequence_2.fa contains the reverse complement from the primer to almost the extremity of the desired family
  ##Packages needed: bbmap and cutadapt, can be both installed by bioconda
  ##First part keeps only the reads where the sequence is present and trim the adapter and part of the transposon
  
  
echo "obtaining informative reads and PCR duplicate removal"
cutadapt --discard-untrimmed -g file:$scriptDir/${te_family}_flanking_sequence.fa -e 0.2 -o $outTrim/sub_1.fastq -p $outTrim/sub_2.fastq $read_1 $read_2 -O 25 -j 8  > $outTrim/trimming_1.txt

  ##This second part remove the reads aligning to the transposon from the P7 read to improve the mapping
cutadapt -A file:$scriptDir/${te_family}_flanking_sequence_2.fa  -e 0.1 -o $outTrim/sub_final_1.fastq -p $outTrim/sub_final_2.fastq $outTrim/sub_1.fastq $outTrim/sub_2.fastq -j 6  > $outTrim/trimming_2.txt
##This third part remove the duplicated reads, based on the P7 read. 
clumpify.sh in=$outTrim/sub_final_2.fastq out=$outTrim/sub_dedup_2.fastq subs=2 dedupe -Xmx32g overwrite=t 2>$outTrim/dedup_1.txt
##Reconstruction of the Paired-end reads used for the analysis
repair.sh in=$outTrim/sub_final_1.fastq in2=$outTrim/sub_dedup_2.fastq out=$outTrim/final_dedup_1.fq.gz out2=$outTrim/final_dedup_2.fq.gz outs=$outTrim/duplicates.fq.gz overwrite=t -Xmx32g 2> $outTrim/dedup_2.txt
echo "informative reads obtained and deduplicated"

rm $outTrim/sub_1.fastq
rm $outTrim/sub_2.fastq
rm $outTrim/sub_final_1.fastq
rm $outTrim/sub_final_2.fastq
rm $outTrim/duplicates.fq.gz
rm $outTrim/sub_dedup_2.fastq

##### MAPPING reads to the reference genome ##### 
bowtie2 -x $refDir -1 $outTrim/final_dedup_1.fq.gz -2 $outTrim/final_dedup_2.fq.gz  -S $outDir/$file/${file}_clip_disc-local.sam --local --very-sensitive -p $CORES --quiet
echo "mapping done"
echo "TE insertion detection"
#Discard reads with a mapping quality lower than 3, can be lowered if you want to recover insertions in highly repetitive sequences alhtough it may increase the number of false positive insertions
samtools view -h -q 3 --threads $CORES $outDir/$file/${file}_clip_disc-local.sam | samtools sort  --threads $CORES -  > $outDir/$file/${file}_clip_disc-local.bam
samtools index  $outDir/$file/${file}_clip_disc-local.bam
rm  $outDir/$file/${file}_clip_disc-local.sam

##### EXTRACT START STOP AND COVERAGE AT ANY POSITION OF THE GENOME ##### 

bamToBed -i  $outDir/$file/${file}_clip_disc-local.bam -split -bed12 | awk ' $1 !="ChrM" && $1 != "ChrC" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$11}' > $outDir/$file/${file}_clip_disc-local_bamtobed.bed


sort -k1,1 -k4 $outDir/$file/${file}_clip_disc-local_bamtobed.bed > $outDir/$file/${file}_clip_disc-local_bamtobed_sorted.bed 

rm $outDir/$file/${file}_clip_disc-local_bamtobed.bed

awk 'BEGIN {
    interval_start=$2 ;interval_stop=$3 ; chr=$1;  name=substr($4,1,index($4,"/")-1) ; qscore=$5; signe=$6; size=$7
    }
    substr($4,1,index($4,"/")-1)==name && $6=="-" {interval_stop=$3; signe=$6;size=$7 ; next}
    substr($4,1,index($4,"/")-1)==name && $6=="+" {interval_start=$2; signe=$6;size=$7 ; next}
    {print chr,interval_start,interval_stop,name,qscore,signe,size,$6,$7 ;  interval_start=$2 ;interval_stop=$3 ; chr=$1;  name=substr($4,1,index($4,"/")-1) ; qscore=$5; signe=$6; size=$7 }' $outDir/$file/${file}_clip_disc-local_bamtobed_sorted.bed   > $outDir/$file/${file}_temp

awk ' {if ($9>=$7 ) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$9 } else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }}' $outDir/$file/${file}_temp | awk 'NR>1  && $3-$2>0 && $3-$2<600 ' | sort -k1,1 -k2,2n >  $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bed


bedToBam  -i  $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bed -g /mnt/data2/leduque/ONT/reference/chrom.sizes > $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bam
samtools index $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bam
samtools depth $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bam -d 5000000 >  $outDir/$file/${file}_depth.bed 
awk '{{if ($6=="-")  print $1"\t"$2+1"\t"$6} if( $6=="+")  print $1"\t"$3"\t"$6}' $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bed | awk '{count[$1" "$2" "$3]++} END {for (word in count) print word, count[word]}' | sort -k1,1 -k2,2n > $outDir/$file/${file}_stop_site.bed
awk '{{if ($6=="-")  print $1"\t"$3"\t"$6} if( $6=="+")  print $1"\t"$2+1"\t"$6}' $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bed | awk '{count[$1" "$2" "$3]++} END {for (word in count) print word, count[word]}' | sort -k1,1 -k2,2n > $outDir/$file/${file}_start_site.bed

#R script
# calculate the ratio (nbr_start+nbr_stop)/depth for all position
Rscript $scriptDir/ratio.R $file $outDir/$file


# ratio 
# "Chr","start","depth","strand","nbr_start"


awk '{if ($5>=$7) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8} if ($5<=$7) { print $1"\t"$2"\t"$3"\t"$6"\t"$5"\t"$7"\t"$8}} ' $outDir/$file/${file}_ratio.bed > $outDir/$file/${file}_TEDseq_ratio_filtered_sorted.bed



rm $outDir/$file/${file}_ratio.bed

##### FILTERING PUTATTIVE INSERTION ##### 




##### DONOR COPY COVERAGE ####### 
### Get the coverage on donor copy to trick the MINCOV parameter:


#intersect with targeted_TE_sequences.bed  to remove parental copy
# Insertion with many many reads ends up in two or three insertion, here we merge them keeping only the one containing the bigest number of read 

awk    -v MINRATIO=$min_ratio -v MINCOV=$min_cov   '$3>=MINCOV &&  $7>MINRATIO '  $outDir/$file/${file}_TEDseq_ratio_filtered_sorted.bed >  $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_sincov.bed

while read line 
 do
 echo $line  | awk '{{if ($4=="+") print $1"\t"$2-400"\t"$2+10"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7} if ($4=="-")  print $1"\t"$2-10"\t"$2+400"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' > $outDir/$file/temp.bed
 awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4}' $outDir/$file/${file}_start_site.bed |  bedtools intersect -a $outDir/$file/temp.bed -b -  -wa -wb | awk '{ sum += $13} END { print sum}' >> $outDir/$file/temp_cov.txt
done <  $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_sincov.bed

paste $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_sincov.bed $outDir/$file/temp_cov.txt > $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_cov.bed
rm  $outDir/$file/temp.bed
rm $outDir/$file/temp_cov.txt

awk  ' {print $1"\t"$2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_cov.bed | bedtools intersect  -v  -a -  -b $ref_TE >$outDir/$file/${file}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed


awk 'BEGIN {
    chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9  }
   $1==chromosome  && $2<=interval_start+40 && $4 < depth { next}
    {print chromosome"\t"interval_start"\t"interval_end"\t"depth"\t"signe"\t"nbr_start"\t"nbr_stop"\t"ratio"\t"cov_pic  ; chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9}'  $outDir/$file/${file}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed  > $outDir/$file/${file}_temp

tail -1 $outDir/$file/${file}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed  >> $outDir/$file/${file}_temp
sort -k1,1 -k2,2nr  $outDir/$file/${file}_temp | tail -1 > $outDir/$file/${file}_temp_2nd
sort -k1,1 -k2,2nr $outDir/$file/${file}_temp | awk 'BEGIN {
    chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9}
    $1==chromosome  && $2>=interval_start-40 && $4 < depth { next}
    {print chromosome"\t"interval_start"\t"interval_end"\t"depth"\t"signe"\t"nbr_start"\t"nbr_stop"\t"ratio"\t"cov_pic  ; chromosome=$1 ; interval_start=$2; interval_end=$3; depth=$4; signe=$5; nbr_start=$6; nbr_stop=$7; ratio=$8; cov_pic=$9}'  >> $outDir/$file/${file}_temp_2nd
sort -k1,1 -k2.2n $outDir/$file/${file}_temp_2nd  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'  |  sed '/^\s*$/d' >   $outDir/$file/${file}_TEDseq_germline_insertion.bed

            rm $outDir/$file/${file}_temp
            rm $outDir/$file/${file}_temp_2nd
            rm $outDir/$file/${file}_TEDseq_epiRILs_insertion_based_on_ratio_filtered.bed 
            rm $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_cov.bed
            rm $outDir/$file/${file}_TEDseq_ratio_filtered_sorted_sincov.bed
            rm $outDir/$file/${file}_depth.bed
            rm $outDir/$file/${file}_start_site.bed
            rm $outDir/$file/${file}_stop_site.bed
            rm $outDir/$file/${file}_TEDseq_ratio_filtered_sorted.bed
            rm $outDir/$file/${file}_clip_disc-local_bamtobed_sorted.bed

mergeBed  -c 4,6,2,3 -o count_distinct,distinct,count_distinct,count_distinct -i  $outDir/$file/${file}_clip_disc-local_bamtobed_nomate.bed | awk '$4>=2 && $3-$2<700' > $outDir/$file/${file}_TED-seq_putative_insertions.bed
awk -v LIBSIZE=$libsize '   $3-$2>LIBSIZE  {print $0}'   $outDir/$file/${file}_TED-seq_putative_insertions.bed | sort -k1,1 -k2,2n | awk -v MINCOV=$min_cov    ' $4>=MINCOV && $6>5 || $7>5  {print $0}' > $outDir/$file/${file}_clean_putative_ins.bed
rm $outDir/$file/${file}_TED-seq_putative_insertions.bed

bedtools intersect -a  $outDir/$file/${file}_TEDseq_germline_insertion.bed  -b  $outDir/$file/${file}_clean_putative_ins.bed  -wa -wb | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$13}' >  $outDir/$file/${file}_TEDseq_germline_insertion_corrected.bed
