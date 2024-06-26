version 1.0

import "imports/pull_smallBamQc_1_0_0.wdl" as smallBamQc

struct FastqInput {
    File read1
    File read2
    String readGroup
}

struct GenomeResources {
    String indexModule 
    String index
    String genomeModule
    String fasta
    String bed
}

workflow emSeqQc {
    input {
        Array[FastqInput] fastqInput
        Int opticalDuplicatePixelDistance
        String outputFileNamePrefix
        String reference
    }

    Map[String, GenomeResources] resources = {
        "hg38": {
            "indexModule": "hg38-bwa-meth-index/p12-2022-10-17",
            "index": "$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa",
            "genomeModule": "hg38-em-seq/p12-2022-10-17",
            "bed": "$HG38_EM_SEQ_ROOT/hg38_random.bed",
            "fasta": "$HG38_EM_SEQ_ROOT/hg38_random.fa"
        }
    }

    GenomeResources ref = resources[reference]

    parameter_meta {
        fastqInput: "A list of Read1 and Read2 FastQs and their readgroup"
        opticalDuplicatePixelDistance: "For MarkDuplicates. The maximum offset between two duplicate clusters in order to consider them optical duplicates. 100 is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate."
        outputFileNamePrefix: "File prefix"
        reference: "Which reference to align to"
    }

    scatter (fastq in fastqInput) {
        call trim_and_align {
            input:
                read1 = fastq.read1,
                read2 = fastq.read2,
                bwaReadGroup = fastq.readGroup,
                bwaIndex = ref.index,
                modules = "fastp/0.23.2 bwa-meth/0.2.5 ~{ref.indexModule}"
        }
    }

    if (length(trim_and_align.bam) == 1) {
        File first_bam = trim_and_align.bam[0]
        File first_index = trim_and_align.index[0]
    }

    # Merging bams won't be done if there is only one bam
    if (length(trim_and_align.bam) > 1) {
        call mergeBams {
            input:
                bams = trim_and_align.bam
        }
    }

    File input_bam = select_first([first_bam, mergeBams.bam])
    File input_bam_index = select_first([first_index, mergeBams.index])

    call mergeFastpJson {
        input:
            jsons = trim_and_align.fastpReport,
            prefix = outputFileNamePrefix
    }
    
    call methylDackel {
        input:
            bam = input_bam,
            index = input_bam_index,
            prefix = outputFileNamePrefix,
            fasta = ref.fasta,
            modules = "methyldackel/0.6.1 ~{ref.genomeModule}"
    }

    call smallBamQc.smallBamQc {
        input:
            bam = input_bam,
            opticalDuplicatePixelDistance = opticalDuplicatePixelDistance,
            outputFileNamePrefix = outputFileNamePrefix
    }

    call samtoolsStatsLambdaControl {
        input:
            bam = input_bam,
            index = input_bam_index,
            prefix = outputFileNamePrefix
    }

    call samtoolsStatsPuc19Control {
        input:
            bam = input_bam,
            index = input_bam_index,
            prefix = outputFileNamePrefix
    }

    output {
        File bedgraph = methylDackel.out
        File fastpReport = mergeFastpJson.json
        File samtools = smallBamQc.samtools
        File picard = smallBamQc.picard
        File bedtoolsCoverage = smallBamQc.bedtoolsCoverage
        File downsampledCounts = smallBamQc.downsampledCounts
        File controlstatsLambda = samtoolsStatsLambdaControl.out
        File controlstatsPuc19 = samtoolsStatsPuc19Control.out
    }


    meta {
        author: "Savo Lazic"
        email: "savo.lazic@oicr.on.ca"
        description: "Generates QC metrics from EM-Seq FastQ input.\n\nAccepts multiple FastQ files, which are trimmed and aligned in parallel. pUC19 and lambda (spiked controls) specific metric files are generated."
        dependencies: [
                      {
                        name: "fastp/0.23.2",
                        url: "https://github.com/OpenGene/fastp"
                      },
                      {
                        name: "bwa-meth/0.2.5",
                        url: "https://github.com/brentp/bwa-meth"
                      },
                      {
                        name: "hg38-bwa-meth-index/p12-2022-10-17",
                        url: "https://gitlab.oicr.on.ca/ResearchIT/modulator"
                      },
                      {
                        name: "methyldackel/0.6.1",
                        url: "https://github.com/dpryan79/MethylDackel"
                      },
                      {
                        name: "samtools/1.15",
                        url: "http://www.htslib.org/doc/samtools.html"
                      }
                      ]
    output_meta: {
    bedgraph: {
        description: "MethylDackel zipped output",
        vidarr_label: "bedgraph"
    },
    fastpReport: {
        description: "Merged fastp json reports",
        vidarr_label: "fastpReport"
    },
    samtools: {
        description: "Samtools stats output",
        vidarr_label: "samtools"
    },
    picard: {
        description: "Picard MarkDuplicates output",
        vidarr_label: "picard"
    },
    bedtoolsCoverage: {
        description: "Bedtools coverage histogram output",
        vidarr_label: "bedtoolsCoverage"
    },
    downsampledCounts: {
        description: "JSON file recording what downsampling was done",
        vidarr_label: "downsampledCounts"
    },
    controlstatsLambda: {
        description: "samtools stats for lambda control only",
        vidarr_label: "controlstatsLambda"
    },
    controlstatsPuc19: {
        description: "samtools stats for pUC19 control only",
        vidarr_label: "controlstatsPuc19"
    }
}
    }
}

task trim_and_align {
    input {
        File read1
        File read2

        String bwaReadGroup
        String bwaIndex

        Boolean fastpDisableQualityFiltering = false
        Int? fastpQualifiedQualityPhred
        Int? fastpUnqualifiedPercentLimit
        Int? fastpNBaseLimit

        Boolean fastpDisableLengthFiltering = false
        Int? fastpLengthRequired

        Boolean fastpDisableAdapterTrimming = false

        Boolean fastpDisableTrimPolyG = false

        Int timeout = 48
        Int memory = 32
        Int threads = 8
        String modules
    }

    parameter_meta {
        read1: "Read 1 FastQ file"
        read2: "Read 2 FastQ file"
        bwaReadGroup: "Read group that will populate the `@RG` BAM flag"
        bwaIndex: "The FastA in the directory that contains the bwa index files"
        fastpDisableQualityFiltering: "Disable fastp quality filtering"
        fastpQualifiedQualityPhred: "The quality value that a base is considered qualified (default >=Q15)"
        fastpUnqualifiedPercentLimit: "How many percents of bases are allowed to be unqualified (default 40%)"
        fastpNBaseLimit: "How many N can a read have before being discarded (default 5)"
        fastpDisableLengthFiltering: "Disable filtering reads below a certain length"
        fastpLengthRequired: "Reads shorter than length_required will be discarded (default 15)"
        fastpDisableAdapterTrimming: "Disable all adapter trimming"
        fastpDisableTrimPolyG: "Disable triming polyG at the end of the read"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    String fastpQ = if fastpDisableQualityFiltering then "-Q" else ""
    String fastpq = if defined(fastpQualifiedQualityPhred) then "-q ~{fastpQualifiedQualityPhred}" else ""
    String fastpu = if defined(fastpUnqualifiedPercentLimit) then "-u ~{fastpUnqualifiedPercentLimit}" else ""
    String fastpn = if defined(fastpNBaseLimit) then "-n ~{fastpNBaseLimit}" else ""

    String fastpL = if fastpDisableLengthFiltering then "-L" else ""
    String fastpl = if defined(fastpNBaseLimit) then "-l ~{fastpNBaseLimit}" else ""

    String fastpA = if fastpDisableAdapterTrimming then "-A" else ""

    String fastpG = if fastpDisableTrimPolyG then "-G" else ""

    command <<<
        set -euo pipefail
        fastp \
            --stdout --thread ~{threads} \
            ~{fastpQ} ~{fastpq} ~{fastpu} ~{fastpn} ~{fastpL} ~{fastpl} ~{fastpA} ~{fastpG} \
            -i ~{read1} -I ~{read2} \
        | bwameth.py -p --threads ~{threads} --read-group ~{bwaReadGroup} --reference ~{bwaIndex} /dev/stdin \
        | samtools sort -o output.bam -@ ~{threads} -
        samtools index -@ ~{threads} output.bam
    >>>

    output {
        File fastpReport = "fastp.json"
        File bam = "output.bam"
        File index = "output.bam.bai"
    }

    meta {
        output_meta: {
            fastpReport: "The json report file produced by fastp",
            bam: "The bam file produced by the trimmed FastQ files fed to bwa-meth",
            index: "The indexed .bai file for the bam"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task mergeBams {
    input {
        Array[File] bams

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    parameter_meta {
        bams: "The bam files to merge"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        set -euo pipefail
        samtools merge -c -p -o output.bam -@ ~{threads} ~{sep=" " bams}
        samtools index -@ ~{threads} output.bam
    >>>

    output {
        File bam = "output.bam"
        File index = "output.bam.bai"
    }

    meta {
        output_meta: {
            bam: "The merged bam",
            index: "The indexed .bai file for the merged bam"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task mergeFastpJson {
    input {
        Array[File] jsons
        String prefix

        Int timeout = 2
        Int memory = 2
        Int threads = 1
        String modules = "jq/1.6"
    }

    parameter_meta {
        jsons: "The fastp json files that will be merged"
        prefix: "File prefix"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        jq -n '[inputs]' ~{sep=" " jsons} > ~{prefix}.fastp.json
    >>>

    output {
        File json = "~{prefix}.fastp.json"
    }

    meta {
        output_meta: {
            json: "The merged json report files"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task methylDackel {
    input {
        File bam
        File index
        String prefix
        String fasta

        Int timeout = 6
        Int memory = 8
        Int threads = 8
        String modules
    }

    parameter_meta {
        bam: "The bam file to analyze"
        index: "The .bai index of the bam file"
        prefix: "File prefix"
        fasta: "FastA file used for alignment"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        set -euo pipefail
        MethylDackel extract --mergeContext -@ ~{threads} ~{fasta} ~{bam} -o ~{prefix}.methyldackel
        gzip ~{prefix}.methyldackel_CpG.bedGraph
    >>>

    output {
        File out = "~{prefix}.methyldackel_CpG.bedGraph.gz"
    }

    meta {
        output_meta: {
            out: "The compressed MethylDackel result bedGraph"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

 task samtoolsStatsLambdaControl {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    parameter_meta {
        bam: "The bam file to analyze"
        index: "The .bai index of the bam file"
        prefix: "File prefix"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        set -euo pipefail
        samtools stats -F 256 -@ ~{threads} ~{bam} lambda > ~{prefix}.lambda.controlstats
    >>>

    output {
        File out = "~{prefix}.lambda.controlstats"
    }

    meta {
        output_meta: {
            out: "Samtools stats output for the lambda control"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

 task samtoolsStatsPuc19Control {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    parameter_meta {
        bam: "The bam file to analyze"
        index: "The .bai index of the bam file"
        prefix: "File prefix"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }

    command <<<
        set -euo pipefail
        samtools stats -F 256 -@ ~{threads} ~{bam} pUC19 > ~{prefix}.puc19.controlstats
    >>>

    output {
        File out = "~{prefix}.puc19.controlstats"
    }

    meta {
        output_meta: {
            out: "Samtools stats output for the pUC19 control"
        }
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}
