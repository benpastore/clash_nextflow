process PREPARE_FOR_RNAUP {

    label 'low'

    publishDir "$params.results/RNAup_prepare/", mode : 'copy', pattern : "*.tsv"
       
    input :
        tuple val(sampleID), val(tsv)

    output :
        tuple val(sampleID), path(".RNAup.tsv"), emit : merge_tsv


    script :
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${tsv} .tsv)
    python3 ${params.bin}/prepare_for_rnaup.py -bwa ${bwa} -tailor ${tailor} -output ${name}.RNAup.tsv

    
    """
}