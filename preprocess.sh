#!/bin/bash
set -x
picardPATH=/data/home/zxy/biosoft
humanGenomeDir=/data/home/zxy/human_genome
workDIR=/data/home/zxy/data/workdir
tmpDIR=/data/home/zxy/data/tmpdir
bwaPATH=/data/home/zxy/miniconda3/bin
dbsnp=/data/home/share/dbBack/knownSNPdb/dbsnp_146.hg38.vcf
oneThousand=/data/home/share/dbBack/knownSNPdb/1000G_phase1.snps.high_confidence.hg38.vcf
millsand=/data/home/share/dbBack/knownSNPdb/Mills_and_1000G_gold_standard.indels.hg38.vcf
GATK_PATH=/data/home/zxy/gatk-4.1.3.0
interval=/data/home/zxy/data/callcnv_GATK/S07604514_Covered_interval.bed

##3.1convert file to sam using picard, all the samples are from the same run, so the read group name is the same
fastq2sam_fun(){
	java -Xmx8G -jar $picardPATH/picard.jar FastqToSam \
		FASTQ=$1 \
		FASTQ2=$2 \
		OUTPUT=$3 \
		READ_GROUP_NAME=A \
		SAMPLE_NAME=$4 \
		LIBRARY_NAME=$4 \
		PLATFORM=illumina 
}

##3.2mark adapter sequences using MarkilluminaAdapters

MarkIlluminaAdapters_fun(){
	java -Xmx8G -jar $picardPATH/picard.jar MarkIlluminaAdapters \
		INPUT=$1 \
		OUTPUT=./${1%_fastqtosam*}_markilluminaadapters.bam \
		METRICS=${1%_fastqtosam*}_markilluminaadapters_metrics.txt \
		TMP_DIR=$workDIR/result/uSamFile
}

##3.3 pipline 

pip_function(){
	java -Xmx8G -jar $picardPATH/picard.jar SamToFastq \
		INPUT=$1 \
		FASTQ=$tmpDIR/stdout \
		CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
		TMP_DIR=$workDIR/tmp | \
		$bwaPATH/bwa mem -M -t 10 -p $humanGenomeDir/Homo_sapiens_assembly38.fasta $tmpDIR/stdin | \
		java -Xmx16G -jar $picardPATH/picard.jar MergeBamAlignment \
		ALIGNED_BAM=$tmpDIR/stdin \
		UNMAPPED_BAM=${1%_markillumina*}_fastqtosam.bam \
		OUTPUT=${1%_markillumina*}_piped.bam \
		R=$humanGenomeDir/Homo_sapiens_assembly38.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
		TMP_DIR=$workDIR/tmp
}

##(3.3.2) change RG id

AddOrReplaceReadGroups_fun(){
	name=${1%.bam*}
	java -Xmx8G -jar $picardPATH/picard.jar AddOrReplaceReadGroups \
		INPUT=$1 \
		OUTPUT=${name}_RG.bam \
		RGID=lane$2 \
		RGLB=${name##*/} \
		RGPL=illumina \
		RGPU=unit1 \
		RGSM=${name##*/} 
	}

#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ001*piped.bam 2
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ002*piped.bam 2
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ003*piped.bam 2
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ004*piped.bam 2
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ005*piped.bam 4
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ006*piped.bam 4
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ007*piped.bam 4
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ008*piped.bam 2
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ009*piped.bam 3
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ010*piped.bam 4
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ011*piped.bam 3
#AddOrReplaceReadGroups_fun $workDIR/result/uSamFile/CJ012*piped.bam 3

##3.4 summerize the outline

CollectAlignmentSummaryMetrics_fun(){
	java -jar $picardPATH/picard.jar CollectAlignmentSummaryMetrics \
		R=$humanGenomeDir/Homo_sapiens_assembly38.fasta \
		I=$1 \
		O=$1.AlignSum.txt \
		-nt 10
}



######split to 1-9 and 10-19 two part
#####no need to do that for you can use seq to get the sequence
#for each in `ls $workDIR/result/uSamFile/CJ00[1-9]*piped_RG.bam`
#do
#	CollectAlignmentSummaryMetrics_fun $each 
#	#echo $each
#done

##3.5 MarkDuplicates

MarkDuplicates_fun(){
	name=${1##*/}
	java -Xmx32G -jar $picardPATH/picard.jar MarkDuplicates\
		INPUT=$1 \
		OUTPUT=${name%.bam*}_markduplicates.bam \
		METRICS_FILE=${name%.bam*}_markduplicates_metrics.txt \
		CREATE_INDEX=true \
		TMP_DIR=$workDIR/tmp

}

#for i in `ls $workDIR/result/uSamFile/*piped_RG.bam`
#do
##mv ${i%_RG*}.bai ${i%.bam*}.bai
#MarkDuplicates_fun $i ### no need for merge the duplicate together in this step
#done


##3.6 base recalibration

BaseRecalibrator_fun(){
	$GATK_PATH BaseRecalibrator \
		--java-options "-Xmx32G" \
		-R $humanGenomeDir/Homo_sapiens_assembly38.fasta \
		-I $1 \
		--known-sites $dbsnp \
		--known-sites $oneThousand \
		--known-sites $millsand \
		-L $interval \
		-O ${1%.bam*}_recal_data.table
}



ApplyBQSR_fun(){
	$GATK_PATH ApplyBQSR \
	--java-options "-Xmx32G" \
	-R $humanGenomeDir/Homo_sapiens_assembly38.fasta \
	-I $1 \
	--bqsr-recal-file ${1%.bam*}_recal_data.table \
	-O ${1%.bam*}_recal.bam
}


MarkIlluminaAdapters_fun $1


