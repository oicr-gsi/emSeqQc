version 1.0

import "imports/pull_bamQC.wdl" as bamQC

struct FastqInput {
    File read1
    File read2
    String readGroup
}

workflow emSeqQc {
    input {
        Array[FastqInput] fastqInput
        String outputFileNamePrefix = "output"
    }

    scatter (fastq in fastqInput) {
        call trim_and_align {
            input:
                read1 = fastq.read1,
                read2 = fastq.read2,
                bwaReadGroup = fastq.readGroup
        }
    }
    call mergeBams {
        input:
            bams = trim_and_align.bam
    }
    call mergeFastpJson {
        input:
            jsons = trim_and_align.fastpReport,
            prefix = outputFileNamePrefix
    }
    call methylDackel {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }
    call bamQC.bamQC {
        input:
            bamFile = mergeBams.bam,
            metadata = {},
            bamQCMetrics_modules = "bam-qc-metrics/0.2.5 hg38-em-seq/p12-2022-10-17",
            bamQCMetrics_refFasta = "$HG38_EM_SEQ_ROOT/hg38_random.fa",
            bamQCMetrics_refSizesBed = "$HG38_EM_SEQ_ROOT/hg38_random.bed",
            bamQCMetrics_workflowVersion = "5.0.2"
    }

    call samtoolsFlagstats {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }

    call samtoolsIdxstats {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }

    call samtoolsStats {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }

    call samtoolsStatsControl {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }

    call bedtoolsGenomeCov {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }

    call gcBias {
        input:
            bam = mergeBams.bam,
            index = mergeBams.index,
            prefix = outputFileNamePrefix
    }

    output {
        File bedgraph = methylDackel.out
        File fastpReport = mergeFastpJson.json
        File bamqc = bamQC.result
        File flagstats = samtoolsFlagstats.out
        File idxstats = samtoolsIdxstats.out
        File stats = samtoolsStats.out
        File controlstats = samtoolsStatsControl.out
        File genomeCov = bedtoolsGenomeCov.out
        File gcBiasOut = gcBias.out
        File gcBiasSummary = gcBias.summary
    }


    meta {
        author: "???"
        email: "???"
        description: "What does the workflow do?"
        dependencies: [
                      {
                          name: "package1/0.1",
                          url: "https://url_to_software"
                      },
                      {
                          name: "package2/0.1",
                          url: "https://url_to_software"
                      }
                      ]
    }
}

task trim_and_align {
    input {
        File read1
        File read2

        String bwaReadGroup
        String bwaReference = "$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa"

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
        String modules = "fastp/0.23.2 bwa-meth/0.2.5 hg38-bwa-meth-index/p12-2022-10-17"
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
        | bwameth.py -p -t ~{threads} --read-group ~{bwaReadGroup} --reference ~{bwaReference} /dev/stdin \
        | samtools sort -o output.bam -@ ~{threads} -
    >>>

    output {
        File fastpReport = "fastp.json"
        File bam = "output.bam"
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

    command <<<
        set -euo pipefail
        samtools merge -c -p -o output.bam -@ ~{threads} ~{sep=" " bams}
        samtools index -@ ~{threads} output.bam
    >>>

    output {
        File bam = "output.bam"
        File index = "output.bam.bai"
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

    command <<<
        jq -n '[inputs]' ~{sep=" " jsons} > ~{prefix}.fastp.json
    >>>

    output {
        File json = "~{prefix}.fastp.json"
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
        String fasta = "$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa"

        Int timeout = 6
        Int memory = 8
        Int threads = 8
        String modules = "methyldackel/0.6.1 hg38-bwa-meth-index/p12-2022-10-17"
    }

    command <<<
        set -euo pipefail
        MethylDackel extract --mergeContext -@ ~{threads} ~{fasta} ~{bam} -o ~{prefix}.methyldackel
        gzip ~{prefix}.methyldackel_CpG.bedGraph
    >>>

    output {
        File out = "~{prefix}.methyldackel_CpG.bedGraph.gz"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task samtoolsFlagstats {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    command <<<
        set -euo pipefail
        samtools flagstats -@ ~{threads} -O tsv ~{bam} > ~{prefix}.flagstats.tsv
    >>>

    output {
        File out = "~{prefix}.flagstats.tsv"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task samtoolsIdxstats {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    command <<<
        set -euo pipefail
        samtools view -h -F 256 -@ ~{threads} ~{bam} \
            | samtools idxstats -@ ~{threads} - > ~{prefix}.idxstats.tsv
    >>>

    output {
        File out = "~{prefix}.idxstats.tsv"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task samtoolsStats {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    command <<<
        set -euo pipefail
        samtools stats -F 256 -@ ~{threads} ~{bam} > ~{prefix}.stats
    >>>

    output {
        File out = "~{prefix}.stats"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

 task samtoolsStatsControl {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 4
        Int threads = 8
        String modules = "samtools/1.15"
    }

    command <<<
        set -euo pipefail
        samtools stats -F 256 -@ ~{threads} ~{bam} lambda pUC19> ~{prefix}.controlstats
    >>>

    output {
        File out = "~{prefix}.controlstats"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task bedtoolsGenomeCov {
    input {
        File bam
        File index
        String prefix

        Int timeout = 6
        Int memory = 8
        Int threads = 8
        String modules = "samtools/1.15 bedtools/2.27"
    }

    command <<<
        set -euo pipefail
        samtools view -u -F 256 ~{bam} \
            | bedtools genomecov -ibam - > ~{prefix}.genomecov.tsv
    >>>

    output {
        File out = "~{prefix}.genomecov.tsv"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task gcBias {
    input {
        File bam
        File index
        String prefix
        String fasta = "$HG38_BWA_METH_INDEX_ROOT/hg38_random.fa"
        Int picardMaxMemMb = 3000

        Int timeout = 6
        Int memory = 6
        Int threads = 1
        String modules = "picard/2.21.2 hg38-bwa-meth-index/p12-2022-10-17 rstats/4.0"
    }

    command <<<
        set -euo pipefail
        samtools view -u -F 256 ~{bam} \
        | java -Xmx~{picardMaxMemMb}M \
            -jar ${PICARD_ROOT}/picard.jar \
            CollectGcBiasMetrics \
            INPUT=/dev/stdin \
            IS_BISULFITE_SEQUENCED=true \
            REFERENCE_SEQUENCE=~{fasta} \
            VALIDATION_STRINGENCY=SILENT \
            OUTPUT=~{prefix}.gcbias.output.txt \
            SUMMARY_OUTPUT=~{prefix}.gcbias.summary.txt \
            CHART_OUTPUT=~{prefix}.gcbias.pdf  # Required, but not part of output
    >>>

    output {
        File out = "~{prefix}.gcbias.output.txt"
        File summary = "~{prefix}.gcbias.summary.txt"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}