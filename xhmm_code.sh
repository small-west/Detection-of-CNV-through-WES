#Using XHMM to detect copy number variation (CNVs) in whole-exome sequencing data
#Author:zxy 
#Date: 2019-11-11
set -x
gatkPATH=/data/home/zxy/biosoft/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
interval=/data/home/zxy/data/interval/S07604514_Covered.bed
reference=/data/home/zxy/human_genome/Homo_sapiens_assembly38.fasta

##write bam.list into 3 groups
##1.run GATK for depth of coverage by gatk3.8
Depthofcoverage(){
	java -Xmx32G -jar $gatkPATH/GenomeAnalysisTK.jar \
		-T DepthOfCoverage \
		-I $1 \
		-L $interval \
		-R $reference \
		-dt BY_SAMPLE -dcov 5000 -l INFO \
		--omitDepthOutputAtEachBase --omitLocusTable \
		--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
		--includeRefNSites --countType COUNT_FRAGMENTS \
		-o ${1}.DATA
}


##2.combine GATK depth of coverage output:
Combine(){
	xhmm --mergeGATKdepths -o ./DATA.RD.txt \
		--GATKdepths group1_CJ1-5.list.DATA.sample_interval_summary \
		--GATKdepths group2_CJ6-11.list.DATA.sample_interval_summary \
		--GATKdepths group3_CJ12-17.list.DATA.sample_interval_summary
}

##3.run GATK to calculate GC content of targets for filter GC content < 0.1 || >0.9 (wd:~/data/interval)
GCcotentbyintervals(){
	java -Xmx3072m -jar $gatkPATH/GenomeAnalysisTK.jar \
		-T GCContentByInterval -L $interval \
		-R $reference \
		-o ./DATA.locus_GC.txt
}
filter(){
	cat ./DATA.locus_GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' >./extreme_gc_targets.txt
}

##4. Filter samples and targets and prepare for normalization
## use the XHMM matrix command to process the read-depth matrix and mean-center the targets in preparation for PCA-based normalization
Preparation(){
	xhmm --matrix -r ./DATA.RD.txt --centerData --centerType target \
		-o ./DATA.filtered_centered.RD.txt \
		--outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
		--outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
		--excludeTargets ./extreme_gc_targets.txt \
		--minTargetSize 10 --maxTargetSize 10000 \
		--minMeanTargetRD 10 --maxMeanTargetRD 500 \
		--minMeanSampleRD 25 --maxMeanSampleRD 200 \
		--maxSdSampleRD 150
}

##5. Run PCA on mean-cemtered data
##determine the strongest independent ways
PCA(){
	xhmm --PCA -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA
}


#6.Normalize mean-centered data using PCA information
##remove the strongest signals
##remove the top principal components and reconstruct a normalized read depth matrix
Normalization(){
	xhmm --normalize -r ./DATA.filtered_centered.RD.txt \
		--PCAfiles ./DATA.RD_PCA \
		--normalizeOutput ./DATA.PCA_normalized.txt \
		--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
}

##7. Filters and z-score centers (by sample) the PCA-normalized data:
Zscore(){
	xhmm --matrix -r ./DATA.PCA_normalized.txt \
		--centerData --centerType sample --zScoreData \
		-o ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
		--outputExcludedTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
		--outputExcludedSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
		--maxSdTargetRD 30
}

##8. Filters original read-depth data to be the same as filtered, normalized data:
Filter(){
	xhmm --matrix -r ./DATA.RD.txt \
		--excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
		--excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
		--excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
		--excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
		-o ./DATA.same_filtered.RD.txt
}

##9. Run the HMM Viterbi algorithm to call CNVs in each sample:
HMM(){
	xhmm --discover -p params.txt \
		-r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
		-R DATA.same_filtered.RD.txt \
		-c DATA.xcnv -a DATA.aux_xcnv -s DATA
}

##10. genotype
Genotype(){
	xhmm --genotype -p params.txt \
		-r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
		-R DATA.same_filtered.RD.txt \
		-g DATA.xcnv -F $reference \
		-v DATA.vcf
}
