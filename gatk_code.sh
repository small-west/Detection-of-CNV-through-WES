##Using gatk to detect germline CNVs in WES data via cohort mode
##Author: zxy 
##Date: 2019-11-12

set -x
referencePATH=/data/home/zxy/human_genome

##1.1 Collect raw counts data with PreprocessIntervals and CollectReadCounts
Preprocessintervals(){
	gatk PreprocessIntervals \
		-R $referencePATH/Homo_sapiens_assembly38.fasta \
		-L S07604514_Covered.bed \
		--bin-length 0 \
		-imr OVERLAPPING_ONLY \
		-O targets.preprocessed.interval_list
}

##1.2 Count reads per bin using CollectReadCounts
Collectreadcounts(){
	gatk CollectReadCounts \
		-L ../interval/targets.preprocessed.interval_list \
		-R $referencePATH/Homo_sapiens_assembly38.fasta \
		-imr OVERLAPPING_ONLY \
		-I $1 \
		--format TSV \
		-O ${1%_piped_markduplicates_recal.bam*}.tsv
}


##2.1 AnnotateIntervals with GC content
AnnotateIntervals(){
	gatk AnnotateIntervals \
		-L ../interval/targets.preprocessed.interval_list \
		-R $referencePATH/Homo_sapiens_assembly38.fasta \
		-imr OVERLAPPING_ONLY \
		-O target.annotated.tsv
}


##2.2FilterIntervals based on GC-content and cohort extreme counts
FilterIntervals(){
	gatk FilterIntervals \
		-L ../interval/targets.preprocessed.interval_list \
		--annotated-intervals ../interval/target.annotated.tsv \
		-I cvg/CJ001_NDHE08576.tsv -I cvg/CJ002_NDHE08577.tsv \
		-I cvg/CJ003_NDHE08578.tsv -I cvg/CJ004_NDHE08579.tsv \
		-I cvg/CJ006_NDHE08590.tsv -I cvg/CJ005_NDHE08589.tsv \
		-I cvg/CJ007_NDHE08591.tsv -I cvg/CJ008_NDHE08580.tsv \
		-I cvg/CJ009_NDHE08581.tsv -I cvg/CJ010_NDHE08592.tsv \
		-I cvg/CJ011_NDHE08582.tsv -I cvg/CJ012_NDHE08583.tsv \
		-I cvg/CJ013_NDHE08584.tsv -I cvg/CJ014_NDHE08585.tsv \
		-I cvg/CJ015_NDHE08586.tsv -I cvg/CJ016_NDHE08593.tsv \
		-I cvg/CJ017_NDHE08587.tsv \
		-imr OVERLAPPING_ONLY \
		-O target.cohort.gc.filtered.interval_list
}


##3.1DetermineGermlineContigPloidy in COHORT MODE
DetermineGermlineContigPloidy(){
	gatk DetermineGermlineContigPloidy \
		-L target.cohort.gc.filtered.interval_list \
		--interval-merging-rule OVERLAPPING_ONLY \
		-I cvg/CJ001_NDHE08576.tsv -I cvg/CJ002_NDHE08577.tsv \
		-I cvg/CJ003_NDHE08578.tsv -I cvg/CJ004_NDHE08579.tsv \
		-I cvg/CJ006_NDHE08590.tsv -I cvg/CJ005_NDHE08589.tsv \
		-I cvg/CJ007_NDHE08591.tsv -I cvg/CJ008_NDHE08580.tsv \
		-I cvg/CJ009_NDHE08581.tsv -I cvg/CJ010_NDHE08592.tsv \
		-I cvg/CJ011_NDHE08582.tsv -I cvg/CJ012_NDHE08583.tsv \
		-I cvg/CJ013_NDHE08584.tsv -I cvg/CJ014_NDHE08585.tsv \
		-I cvg/CJ015_NDHE08586.tsv -I cvg/CJ016_NDHE08593.tsv \
		-I cvg/CJ017_NDHE08587.tsv \
		--contig-ploidy-priors target_contig_ploidy_priors.tsv \
		--output . \
		--output-prefix ploidy \
		--verbosity DEBUG
}

##4.1 Split gc.filtered.interval_list into 11 scatters, 20000 bins per scatter. Just metain 10Mbp-50Mbp per scatter.
SplitIntervals(){
	gatk IntervalListTools \
		--INPUT target.cohort.gc.filtered.interval_list \
		--SUBDIVISION_MODE INTERVAL_COUNT \
		--SCATTER_CONTENT 20000 \
		--OUTPUT scatter
}

##4.2 GermlineCNVcaller in cohort mode $1:output of 4.1 $2:code of scatter
GermlineCNVCaller(){
	gatk GermlineCNVCaller \
		--run-mode COHORT \
		-L $1 \
		-I cvg/CJ001_NDHE08576.tsv -I cvg/CJ002_NDHE08577.tsv \
		-I cvg/CJ003_NDHE08578.tsv -I cvg/CJ004_NDHE08579.tsv \
		-I cvg/CJ006_NDHE08590.tsv -I cvg/CJ005_NDHE08589.tsv \
		-I cvg/CJ007_NDHE08591.tsv -I cvg/CJ008_NDHE08580.tsv \
		-I cvg/CJ009_NDHE08581.tsv -I cvg/CJ010_NDHE08592.tsv \
		-I cvg/CJ011_NDHE08582.tsv -I cvg/CJ012_NDHE08583.tsv \
		-I cvg/CJ013_NDHE08584.tsv -I cvg/CJ014_NDHE08585.tsv \
		-I cvg/CJ015_NDHE08586.tsv -I cvg/CJ016_NDHE08593.tsv \
		-I cvg/CJ017_NDHE08587.tsv \
		--contig-ploidy-calls ploidy-calls \
		--annotated-intervals ../interval/target.annotated.tsv \
		--interval-merging-rule OVERLAPPING_ONLY \
		--output cohort17-target \
		--output-prefix cohort17-target${2}_of11 \
		--verbosity DEBUG
}

##5. PostprocessGermlineCNVCalls COHORT MODE $1:index of sample and it begin with 0, $2:Sample name
PostprocessGermlineCNVCalls(){
	gatk PostprocessGermlineCNVCalls \
		--model-shard-path cohort17-target/cohort17-target01_of11-model \
		--model-shard-path cohort17-target/cohort17-target02_of11-model \
		--model-shard-path cohort17-target/cohort17-target03_of11-model \
		--model-shard-path cohort17-target/cohort17-target04_of11-model \
		--model-shard-path cohort17-target/cohort17-target05_of11-model \
		--model-shard-path cohort17-target/cohort17-target06_of11-model \
		--model-shard-path cohort17-target/cohort17-target07_of11-model \
		--model-shard-path cohort17-target/cohort17-target08_of11-model \
		--model-shard-path cohort17-target/cohort17-target09_of11-model \
		--model-shard-path cohort17-target/cohort17-target10_of11-model \
		--model-shard-path cohort17-target/cohort17-target11_of11-model \
		--calls-shard-path cohort17-target/cohort17-target01_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target02_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target03_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target04_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target05_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target06_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target07_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target08_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target09_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target10_of11-calls \
		--calls-shard-path cohort17-target/cohort17-target11_of11-calls \
		--allosomal-contig chrX --allosomal-contig chrY \
		--contig-ploidy-calls ploidy-calls \
		--sample-index $1 \
		--output-genotyped-intervals genotyped-intervals-${2}.vcf.gz \
		--output-genotyped-segments genotyped-segments-${2}.vcf.gz \
		--output-denoised-copy-ratios denoised-copy-ratios-${2}.tsv \
		--sequence-dictionary $referencePATH/Homo_sapiens_assembly38.dict
}






