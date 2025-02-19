process MERGE_TAILOR_BWA {

    label 'low'

    publishDir "$params.results/merge_tailor_bwa_alignment/", mode : 'copy', pattern : "*.tsv"
       
    input :
        tuple val(sampleID), val(bwa)
        tuple val(sampleID), val(tailor)

    output :
        tuple val(sampleID), path("*.merge.tsv"), emit : merge_tsv


    script :
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${tailor} .bed)
    python3 ${params.bin}/merge_tailor_bwa_alignments.py -bwa ${bwa} -tailor ${tailor} -output ${name}.merge.tsv

    
    """
}