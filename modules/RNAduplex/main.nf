process RNADUPLEX {

    label 'low'

    publishDir "$params.results/RNAduplex/", mode : 'copy', pattern : "*.tsv"
       
    input :
        tuple val(sampleID), val(tsv)
        val transcripts
        val smrnas

    output :
        tuple val(sampleID), path("*RNAduplex.tsv"), emit : tsv


    script :
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${tsv} .tsv)

    python3 ${params.bin}/rnaduplex_v5.py ${tsv} ${transcripts} ${smrnas} \$name.RNAduplex.tsv

    """
}