process TAILOR_INDEX {

    label 'high'
    
    publishDir "$params.index/tailor", mode : 'copy'

    input :
        val genome
        val index
        val name

    output :
        val(index), emit : tailor_index
        val("*chrom_sizes"), emit : chrom_sizes
        path("${name}/*")
        val("*chrLen")
        val("*chrStart")
        val("*NposLen.z")
        val("*_bwt.bwt")
        val("*_seq.bwt")
        val("*_table.bwt")

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    [ ! -d ${name} ] && mkdir -p ${name}

    ${params.bin}/tailor_v11 build -i ${genome} -p ${name}/${name}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}/${name}_chrom_sizes
    rm ${genome}.fai
    """
}

process TAILOR_ALIGN {

    label 'high'

    publishDir "$params.results/tailor/", mode : 'copy', pattern : "*.bed"
       
    input :
        val talor_index
        tuple val(sampleID), val(tailor_input)

    output :
        tuple val(sampleID), path("*.tailor.bed"), emit : tailor_bed


    script :
    """
    #!/bin/bash

    source activate smrnaseq
    
    # bwa output is gene, start, end, read_name, score, strand, aligned_portion, unaligned_portion
    # we want to take unaligned portion and align it to the smRNA reference

    name=\$(basename ${tailor_input} .bed)
    cat ${tailor_input} | cut -f4,8 | awk '(\$2!="")' | awk -F'\\t' '{ print ">"\$1"\\n"\$2}' > uniq_unal

    python3 ${params.bin}/fasta_to_fastq.py uniq_unal > uniq_unal.fq


    # analyze 3' nontemplated nucleotides
    ${params.bin}/tailor_v11 map \\
        -i uniq_unal.fq \\
        -p ${talor_index} \\
        -n ${task.cpus} \\
        2> tailor.genome.log | \\
    tee \$sam | \\
    ${params.bin}/tailor_sam_to_bed | \\
    awk '(\$5 <= 2 && \$6 == "+")' | \\
    awk -v num=\$nTag -F'\\t' -v OFS="\\t" '{
        print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8
    }' | \\
    awk '(\$2<=3)' | \\
    awl '(length(\$7) >= 17)' > \$name.tailor.bed

    """
}