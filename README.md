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
`opticalDuplicatePixelDistance`|Int|For MarkDuplicates. The maximum offset between two duplicate clusters in order to consider them optical duplicates. 100 is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate.
`outputFileNamePrefix`|String|File prefix
`reference`|String|Which reference to align to


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
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
`mergeBams.timeout`|Int|6|The hours until the task is killed
`mergeBams.memory`|Int|4|The GB of memory provided to the task
`mergeBams.threads`|Int|8|The number of threads the task has access to
`mergeBams.modules`|String|"samtools/1.15"|The modules that will be loaded
`mergeFastpJson.timeout`|Int|2|The hours until the task is killed
`mergeFastpJson.memory`|Int|2|The GB of memory provided to the task
`mergeFastpJson.threads`|Int|1|The number of threads the task has access to
`mergeFastpJson.modules`|String|"jq/1.6"|The modules that will be loaded
`methylDackel.timeout`|Int|6|The hours until the task is killed
`methylDackel.memory`|Int|8|The GB of memory provided to the task
`methylDackel.threads`|Int|8|The number of threads the task has access to
`smallBamQc.storeDownsampledCounts_timeout`|Int|1|The hours until the task is killed.
`smallBamQc.storeDownsampledCounts_threads`|Int|1|The number of threads the task has access to.
`smallBamQc.storeDownsampledCounts_memory`|Float|0.1|The GB of memory provided to the task.
`smallBamQc.storeDownsampledCounts_modules`|String|""|The modules that will be loaded.
`smallBamQc.bedtoolsCoverageSmall_timeout`|Int|1|The hours until the task is killed.
`smallBamQc.bedtoolsCoverageSmall_threads`|Int|2|The number of threads the task has access to.
`smallBamQc.bedtoolsCoverageSmall_memory`|Int|4|The GB of memory provided to the task.
`smallBamQc.bedtoolsCoverageSmall_modules`|String|"bedtools/2.27 samtools/1.16.1"|The modules that will be loaded.
`smallBamQc.bedtoolsCoverageFull_timeout`|Int|12|The hours until the task is killed.
`smallBamQc.bedtoolsCoverageFull_threads`|Int|2|The number of threads the task has access to.
`smallBamQc.bedtoolsCoverageFull_memory`|Int|4|The GB of memory provided to the task.
`smallBamQc.bedtoolsCoverageFull_modules`|String|"bedtools/2.27 samtools/1.16.1"|The modules that will be loaded.
`smallBamQc.featuresHead_timeout`|Int|1|The hours until the task is killed.
`smallBamQc.featuresHead_threads`|Int|1|The number of threads the task has access to.
`smallBamQc.featuresHead_memory`|Float|0.1|The GB of memory provided to the task.
`smallBamQc.featuresHead_modules`|String|""|The modules that will be loaded.
`smallBamQc.markDuplicates_timeout`|Int|12|hours before task timeout.
`smallBamQc.markDuplicates_threads`|Int|4|Requested CPU threads.
`smallBamQc.markDuplicates_jobMemory`|Int|16|Memory allocated for this job.
`smallBamQc.markDuplicates_modules`|String|"picard/2.21.2"|required environment modules.
`smallBamQc.markDuplicates_picardMaxMemMb`|Int|6000|Memory requirement in MB for running Picard JAR.
`smallBamQc.samtoolsHead_modules`|String|"samtools/1.16.1"|The modules that will be loaded.
`smallBamQc.samtoolsHead_threads`|Int|4|The number of threads the task has access to.
`smallBamQc.samtoolsHead_memory`|Int|2|The GB of memory provided to the task.
`smallBamQc.samtoolsHead_timeout`|Int|6|The hours until the task is killed.
`smallBamQc.samtoolsStatsFull_modules`|String|"samtools/1.16.1"|The modules that will be loaded.
`smallBamQc.samtoolsStatsFull_threads`|Int|4|The number of threads the task has access to.
`smallBamQc.samtoolsStatsFull_memory`|Int|2|The GB of memory provided to the task.
`smallBamQc.samtoolsStatsFull_timeout`|Int|6|The hours until the task is killed.
`smallBamQc.samtoolsStatsSmall_modules`|String|"samtools/1.16.1"|The modules that will be loaded.
`smallBamQc.samtoolsStatsSmall_threads`|Int|4|The number of threads the task has access to.
`smallBamQc.samtoolsStatsSmall_memory`|Int|2|The GB of memory provided to the task.
`smallBamQc.samtoolsStatsSmall_timeout`|Int|1|The hours until the task is killed.
`smallBamQc.bedtoolsReadsToUse`|Int?|None|If defined, use that many reads from the beginning of the BAM file for bedtools analysis. If not defined, use all BAM reads.
`smallBamQc.features`|String?|None|If defined, bedtools calculates coverage for those features only. If not defined, calculate coverage for all reads (whole genome).
`smallBamQc.featuresToUse`|Int?|None|If defined, use that many features from the beginning of the features file. If not defined, use all features.
`smallBamQc.picardMarkDuplicatesReadsToUse`|Int?|None|If defined, MarkDuplicates uses that many reads from the beginning of the BAM file. If not defined, use all BAM reads. Note that a new BAM file is created if defined, so using a large number will temporarily generate a second large BAM file.
`smallBamQc.samtoolsStatsReadsToUse`|Int?|None|If defined, use that many read from the beginning of the BAM file for samtools stats. If not defined, use all BAM reads.
`samtoolsStatsLambdaControl.timeout`|Int|6|The hours until the task is killed
`samtoolsStatsLambdaControl.memory`|Int|4|The GB of memory provided to the task
`samtoolsStatsLambdaControl.threads`|Int|8|The number of threads the task has access to
`samtoolsStatsLambdaControl.modules`|String|"samtools/1.15"|The modules that will be loaded
`samtoolsStatsPuc19Control.timeout`|Int|6|The hours until the task is killed
`samtoolsStatsPuc19Control.memory`|Int|4|The GB of memory provided to the task
`samtoolsStatsPuc19Control.threads`|Int|8|The number of threads the task has access to
`samtoolsStatsPuc19Control.modules`|String|"samtools/1.15"|The modules that will be loaded


### Outputs

Output | Type | Description | Labels
---|---|---|---
`bedgraph`|File|MethylDackel zipped output|vidarr_label: bedgraph
`fastpReport`|File|Merged fastp json reports|vidarr_label: fastpReport
`samtools`|File|Samtools stats output|vidarr_label: samtools
`picard`|File|Picard MarkDuplicates output|vidarr_label: picard
`bedtoolsCoverage`|File|Bedtools coverage histogram output|vidarr_label: bedtoolsCoverage
`downsampledCounts`|File|JSON file recording what downsampling was done|vidarr_label: downsampledCounts
`controlstatsLambda`|File|samtools stats for lambda control only|vidarr_label: controlstatsLambda
`controlstatsPuc19`|File|samtools stats for pUC19 control only|vidarr_label: controlstatsPuc19


## Commands
 This section lists command(s) run by emSeqQC workflow
 
 * Running emSeqQC
 
 
 ```
         set -euo pipefail
         fastp \
             --stdout --thread ~{threads} \
             ~{fastpQ} ~{fastpq} ~{fastpu} ~{fastpn} ~{fastpL} ~{fastpl} ~{fastpA} ~{fastpG} \
             -i ~{read1} -I ~{read2} \
         | bwameth.py -p --threads ~{threads} --read-group ~{bwaReadGroup} --reference ~{bwaIndex} /dev/stdin \
         | samtools sort -o output.bam -@ ~{threads} -
         samtools index -@ ~{threads} output.bam
 ```
 
 ### Merging and Indexing bams
 
 ```
         set -euo pipefail
         samtools merge -c -p -o output.bam -@ ~{threads} ~{sep=" " bams}
         samtools index -@ ~{threads} output.bam
 ```
 ### Merging json reports
 
 ```
         jq -n '[inputs]' ~{sep=" " jsons} > ~{prefix}.fastp.json
 ```
 ### Running MethylDackel
 
 ```
         set -euo pipefail
         MethylDackel extract --mergeContext -@ ~{threads} ~{fasta} ~{bam} -o ~{prefix}.methyldackel
         gzip ~{prefix}.methyldackel_CpG.bedGraph
 ```
 
 ### Running samtools stats, no unmapped or secondary reads. lambda reference
 
 ```
         set -euo pipefail
         samtools stats -F 256 -@ ~{threads} ~{bam} lambda > ~{prefix}.lambda.controlstats
 ```
 
 ### Running samtools stats, no unmapped or secondary reads. pUC19 reference
 
 ```
         set -euo pipefail
         samtools stats -F 256 -@ ~{threads} ~{bam} pUC19 > ~{prefix}.puc19.controlstats
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
