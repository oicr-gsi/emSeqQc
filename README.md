# emSeqQc

Generates QC metrics from EM-Seq FastQ input.

Accepts multiple FastQ files, which are trimmed and aligned in parallel. pUC19 and lambda (spiked controls) specific metric files are generated.

## Overview

## Dependencies

* [fastp 0.23.2](https://github.com/OpenGene/fastp)
* [bwa-meth 0.2.5](https://github.com/brentp/bwa-meth)
* [hg38-bwa-meth-index p12-2022-10-17](https://gitlab.oicr.on.ca/ResearchIT/modulator)
* [methyldackel 0.6.1](https://github.com/dpryan79/MethylDackel)
* [samtools 1.15](http://www.htslib.org/doc/samtools.html)


## Usage

### Cromwell
```
java -jar cromwell.jar run emSeqQc.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`fastqInput`|Array[FastqInput]|A list of Read1 and Read2 FastQs and their readgroup


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|"output"|File prefix


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`trim_and_align.bwaIndex`|String|"$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa"|The FastA in the directory that contains the bwa index files
`trim_and_align.fastpDisableQualityFiltering`|Boolean|false|Disable fastp quality filtering
`trim_and_align.fastpQualifiedQualityPhred`|Int?|None|The quality value that a base is considered qualified (default >=Q15)
`trim_and_align.fastpUnqualifiedPercentLimit`|Int?|None|How many percents of bases are allowed to be unqualified (default 40%)
`trim_and_align.fastpNBaseLimit`|Int?|None|How many N can a read have before being discarded (default 5)
`trim_and_align.fastpDisableLengthFiltering`|Boolean|false|Disable filtering reads below a certain length
`trim_and_align.fastpLengthRequired`|Int?|None|Reads shorter than length_required will be discarded (default 15)
`trim_and_align.fastpDisableAdapterTrimming`|Boolean|false|Disable all adapter trimming
`trim_and_align.fastpDisableTrimPolyG`|Boolean|false|Disable triming polyG at the end of the read
`trim_and_align.timeout`|Int|48|The hours until the task is killed
`trim_and_align.memory`|Int|32|The GB of memory provided to the task
`trim_and_align.threads`|Int|8|The number of threads the task has access to
`trim_and_align.modules`|String|"fastp/0.23.2 bwa-meth/0.2.5 hg38-bwa-meth-index/p12-2022-10-17"|The modules that will be loaded
`mergeBams.timeout`|Int|6|The hours until the task is killed
`mergeBams.memory`|Int|4|The GB of memory provided to the task
`mergeBams.threads`|Int|8|The number of threads the task has access to
`mergeBams.modules`|String|"samtools/1.15"|The modules that will be loaded
`mergeFastpJson.timeout`|Int|2|The hours until the task is killed
`mergeFastpJson.memory`|Int|2|The GB of memory provided to the task
`mergeFastpJson.threads`|Int|1|The number of threads the task has access to
`mergeFastpJson.modules`|String|"jq/1.6"|The modules that will be loaded
`methylDackel.fasta`|String|"$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa"|FastA file used for alignment
`methylDackel.timeout`|Int|6|The hours until the task is killed
`methylDackel.memory`|Int|8|The GB of memory provided to the task
`methylDackel.threads`|Int|8|The number of threads the task has access to
`methylDackel.modules`|String|"methyldackel/0.6.1 hg38-bwa-meth-index/p12-2022-10-17"|The modules that will be loaded
`bamQC.collateResults_timeout`|Int|1|hours before task timeout
`bamQC.collateResults_threads`|Int|4|Requested CPU threads
`bamQC.collateResults_jobMemory`|Int|8|Memory allocated for this job
`bamQC.collateResults_modules`|String|"python/3.6"|required environment modules
`bamQC.cumulativeDistToHistogram_timeout`|Int|1|hours before task timeout
`bamQC.cumulativeDistToHistogram_threads`|Int|4|Requested CPU threads
`bamQC.cumulativeDistToHistogram_jobMemory`|Int|8|Memory allocated for this job
`bamQC.cumulativeDistToHistogram_modules`|String|"python/3.6"|required environment modules
`bamQC.runMosdepth_timeout`|Int|4|hours before task timeout
`bamQC.runMosdepth_threads`|Int|4|Requested CPU threads
`bamQC.runMosdepth_jobMemory`|Int|16|Memory allocated for this job
`bamQC.runMosdepth_modules`|String|"mosdepth/0.2.9"|required environment modules
`bamQC.bamQCMetrics_timeout`|Int|4|hours before task timeout
`bamQC.bamQCMetrics_threads`|Int|4|Requested CPU threads
`bamQC.bamQCMetrics_jobMemory`|Int|16|Memory allocated for this job
`bamQC.bamQCMetrics_normalInsertMax`|Int|1500|Maximum of expected insert size range
`bamQC.markDuplicates_timeout`|Int|4|hours before task timeout
`bamQC.markDuplicates_threads`|Int|4|Requested CPU threads
`bamQC.markDuplicates_jobMemory`|Int|16|Memory allocated for this job
`bamQC.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules
`bamQC.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR
`bamQC.markDuplicates_opticalDuplicatePixelDistance`|Int|100|Maximum offset between optical duplicate clusters
`bamQC.downsampleRegion_timeout`|Int|4|hours before task timeout
`bamQC.downsampleRegion_threads`|Int|4|Requested CPU threads
`bamQC.downsampleRegion_jobMemory`|Int|16|Memory allocated for this job
`bamQC.downsampleRegion_modules`|String|"samtools/1.9"|required environment modules
`bamQC.downsample_timeout`|Int|4|hours before task timeout
`bamQC.downsample_threads`|Int|4|Requested CPU threads
`bamQC.downsample_jobMemory`|Int|16|Memory allocated for this job
`bamQC.downsample_modules`|String|"samtools/1.9"|required environment modules
`bamQC.downsample_randomSeed`|Int|42|Random seed for pre-downsampling (if any)
`bamQC.downsample_downsampleSuffix`|String|"downsampled.bam"|Suffix for output file
`bamQC.findDownsampleParamsMarkDup_timeout`|Int|4|hours before task timeout
`bamQC.findDownsampleParamsMarkDup_threads`|Int|4|Requested CPU threads
`bamQC.findDownsampleParamsMarkDup_jobMemory`|Int|16|Memory allocated for this job
`bamQC.findDownsampleParamsMarkDup_modules`|String|"python/3.6"|required environment modules
`bamQC.findDownsampleParamsMarkDup_customRegions`|String|""|Custom downsample regions; overrides chromosome and interval parameters
`bamQC.findDownsampleParamsMarkDup_intervalStart`|Int|100000|Start of interval in each chromosome, for very large BAMs
`bamQC.findDownsampleParamsMarkDup_baseInterval`|Int|15000|Base width of interval in each chromosome, for very large BAMs
`bamQC.findDownsampleParamsMarkDup_chromosomes`|Array[String]|["chr12", "chr13", "chrXII", "chrXIII"]|Array of chromosome identifiers for downsampled subset
`bamQC.findDownsampleParamsMarkDup_threshold`|Int|10000000|Minimum number of reads to conduct downsampling
`bamQC.findDownsampleParams_timeout`|Int|4|hours before task timeout
`bamQC.findDownsampleParams_threads`|Int|4|Requested CPU threads
`bamQC.findDownsampleParams_jobMemory`|Int|16|Memory allocated for this job
`bamQC.findDownsampleParams_modules`|String|"python/3.6"|required environment modules
`bamQC.findDownsampleParams_preDSMultiplier`|Float|1.5|Determines target size for pre-downsampled set (if any). Must have (preDSMultiplier) < (minReadsRelative).
`bamQC.findDownsampleParams_precision`|Int|8|Number of decimal places in fraction for pre-downsampling
`bamQC.findDownsampleParams_minReadsRelative`|Int|2|Minimum value of (inputReads)/(targetReads) to allow pre-downsampling
`bamQC.findDownsampleParams_minReadsAbsolute`|Int|10000|Minimum value of targetReads to allow pre-downsampling
`bamQC.findDownsampleParams_targetReads`|Int|100000|Desired number of reads in downsampled output
`bamQC.indexBamFile_timeout`|Int|4|hours before task timeout
`bamQC.indexBamFile_threads`|Int|4|Requested CPU threads
`bamQC.indexBamFile_jobMemory`|Int|16|Memory allocated for this job
`bamQC.indexBamFile_modules`|String|"samtools/1.9"|required environment modules
`bamQC.countInputReads_timeout`|Int|4|hours before task timeout
`bamQC.countInputReads_threads`|Int|4|Requested CPU threads
`bamQC.countInputReads_jobMemory`|Int|16|Memory allocated for this job
`bamQC.countInputReads_modules`|String|"samtools/1.9"|required environment modules
`bamQC.updateMetadata_timeout`|Int|4|hours before task timeout
`bamQC.updateMetadata_threads`|Int|4|Requested CPU threads
`bamQC.updateMetadata_jobMemory`|Int|16|Memory allocated for this job
`bamQC.updateMetadata_modules`|String|"python/3.6"|required environment modules
`bamQC.filter_timeout`|Int|4|hours before task timeout
`bamQC.filter_threads`|Int|4|Requested CPU threads
`bamQC.filter_jobMemory`|Int|16|Memory allocated for this job
`bamQC.filter_modules`|String|"samtools/1.9"|required environment modules
`bamQC.filter_minQuality`|Int|30|Minimum alignment quality to pass filter
`bamQC.outputFileNamePrefix`|String|"bamQC"|Prefix for output files
`samtoolsStatsLambdaControl.timeout`|Int|6|The hours until the task is killed
`samtoolsStatsLambdaControl.memory`|Int|4|The GB of memory provided to the task
`samtoolsStatsLambdaControl.threads`|Int|8|The number of threads the task has access to
`samtoolsStatsLambdaControl.modules`|String|"samtools/1.15"|The modules that will be loaded
`samtoolsStatsPuc19Control.timeout`|Int|6|The hours until the task is killed
`samtoolsStatsPuc19Control.memory`|Int|4|The GB of memory provided to the task
`samtoolsStatsPuc19Control.threads`|Int|8|The number of threads the task has access to
`samtoolsStatsPuc19Control.modules`|String|"samtools/1.15"|The modules that will be loaded


### Outputs

Output | Type | Description
---|---|---
`bedgraph`|File|MethylDackel zipped output
`fastpReport`|File|Merged fastp json reports
`bamqc`|File|bamqc json
`controlstatsLambda`|File|samtools stats for lambda control only
`controlstatsPuc19`|File|samtools stats for pUC19 control only


## Commands
 See WDL.
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
