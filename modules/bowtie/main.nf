process BOWTIE_INDEX_GENOME {

    label 'low'
    
    publishDir "$params.index/bowtie", mode : 'copy'

    input :
        val genome
        val index
        val chrom_sizes
        val name

    output :
        val(index), emit : bowtie_index
        val(chrom_sizes), emit : bowtie_chrom_sizes
        path("${name}/*")

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    [ ! -d ${name} ] && mkdir -p ${name}

    # bowtie
    bowtie-build --quiet ${genome} ${name}/${name} --threads ${task.cpus}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}/${name}_chrom_sizes
    rm ${genome}.fai
    """
}

process BOWTIE_ALIGN_GENOME {

    label 'low'

    publishDir "$params.results/bowtie", mode : 'copy', pattern : '*'

    input : 
        val idx
        tuple val(sampleID), path(fastq)

    output :
        tuple val(sampleID), path("unmapped.fq"), emit : unmapped
        path("bowtie_alignment_stats.txt")

    script : 
    """
    #!/bin/bash

    source activate smrnaseq

    #############################################################################
    # Do Alignment with -a --best --strata, get all aligned reads in best stratum
    #############################################################################
    cat $fastq | \\
        bowtie \\
        -x ${idx} \\
        -q - \\
        -p ${task.cpus} \\
        -a \\
        --un unmapped.fq \\
        --best \\
        --strata \\
        -v 3 \\
        -m 1000 \\
        -S > aligned.sam 2> bowtie_alignment_stats.txt

    """
}

process BOWTIE_ALIGN_CLASH_TARGET_SITE_GENOME {

    label 'low'

    publishDir "$params.results/bowtie", mode : 'copy', pattern : '*.bed'

    input : 
        val idx
        tuple val(sampleID), path(clash_output)

    output :
        tuple val(sampleID), path("*.bed"), emit : unmapped
        path("bowtie_align_clash_target_site_genome.txt")

    script : 
    """
    #!/bin/bash

    source activate smrnaseq



    #############################################################################
    # Do Alignment with -a --best --strata, get all aligned reads in best stratum
    #############################################################################
    cat fasta | \\
        bowtie \\
        -x ${idx} \\
        -q - \\
        -p ${task.cpus} \\
        -a \\
        --un unmapped.fq \\
        --best \\
        --strata \\
        -v 3 \\
        -m 1000 \\
        -S > aligned.sam 2> bowtie_align_clash_target_site_genome.txt

    """
}