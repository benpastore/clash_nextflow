process TRIM_GALORE {

    label 'low'

    publishDir "$params.results/trim_galore/logs", mode : 'copy', pattern : "*_trimming_report.txt"
    publishDir "$params.results/trim_galore/collapsed", mode : 'copy', pattern : "*uniq.fastq"
    publishDir "$params.results/trim_galore/fastqc", mode : 'copy', pattern : '*_fastqc.{zip,html}'

    input :
        tuple val(sampleID), val(fastq)
    
    output : 
        path("${sampleID}*_trimming_report.txt")
        tuple val(sampleID), path("uniq.fastq"), emit : trimmed_fq
        //tuple val(sampleID), path("${sampleID}.trimmed.uniq.fa"), emit : collapsed_fa
        path("*_fastqc.{zip,html}")

    script : 
    adapter = params.adapter ? "-a ${params.adapter}" : ''
    min_length = params.min_length ? "--length ${params.min_length}" : ''
    max_length = params.max_length ? "--max_length ${params.max_length}" : ''
    quality = params.quality || params.quality == 0 ? "-q ${params.quality}" : '-q 30'
    """
    #!/bin/bash 

    source activate smrnaseq

    fq=${fastq}
    t_fq=${sampleID}_trimmed.fq.gz
    fa=${sampleID}.trimmed.uniq.fa

    # run fastqc on data before trimming
    fastqc \$fq -o \$PWD

    # trim adapters
    trim_galore -j ${task.cpus} ${adapter} ${quality} -e 0.1 --gzip ${min_length} ${max_length} --fastqc \$fq

    python3 ${params.bin}/fastq_to_uni_fastq.py \$fq uniq.fastq

    
    """
}