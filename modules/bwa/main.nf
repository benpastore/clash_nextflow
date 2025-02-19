process BWA_INDEX {

    tag "${genome}_BWA_INDEX"

    label 'bwa_low'

    publishDir "$params.index/bwa", mode : 'copy'
    
    input : 
        val genome
        val index
        val name 

    output : 
        val(index), emit : bwa_index
        path("${name}/*")

    script :
    """
    #!/bin/bash

    source activate rnaseq

    bwa index -p ${name} ${genome}
    
    samtools faidx ${genome}
    cut -f1,2 ${genome}.fai > chrom_sizes.txt

    mkdir ${name}
    find . -type f -name "*${name}*" -exec cp {} ${name} \\;
    #cp ${name}* ${name}
    cp chrom_sizes.txt ${name}
    """
}

process BWA_MEM {

    tag "${sampleID}_BWA_MEM"

    label 'bwa_high'

    publishDir "$params.results/bwa/", mode : 'copy', pattern : "*.txt"
    publishDir "$params.results/bwa/", mode : 'copy', pattern : "*.sam"

    input : 
        val bwa_idx
        tuple val(sampleID), val(fastq)

    output : 
        tuple val(sampleID), path("${sampleID}.bwa.sam"), emit : bwa_sam
        path("${sampleID}.txt")
        
    script :
    """
    #!/bin/bash

    source activate rnaseq

    bwa mem ${bwa_idx} \\
        -t ${task.cpus} \\
        -T 10 \\
        ${fastq} > ${sampleID}.bwa.sam

    samtools flagstat ${sampleID}.bwa.sam > ${sampleID}.txt

    """
}

process PROCESS_BWA_MEM {

    label 'low'

    publishDir "$params.results/bwa/", mode : 'copy', pattern : "*.bed"

    input : 
        tuple val(sampleID), val(sam)
    
    output : 
        tuple val(sampleID), path("*processed.bed"), emit: bwa_bed
    
    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    name=\$(basename ${sam} .sam)

    python3 ${params.bin}/process_bwa_output.py ${sam} \$name.processed.bed
    """

}