process BOWTIE2_INDEX_GENOME {

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
    module load bowtie2/2.4.1

    [ ! -d ${name} ] && mkdir -p ${name}

    # bowtie
    bowtie2-build --quiet ${genome} ${name}/${name} --threads ${task.cpus}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}/${name}_chrom_sizes
    rm ${genome}.fai
    """
}

process BOWTIE2_INDEX_JUNCTION {

    label 'low'
    
    publishDir "$params.index/bowtie", mode : 'copy'

    input :
        val juncs
        val index
        val name

    output :
        val(index), emit : bowtie_index
        path("${name}/*")

    script :
    """
    #!/bin/bash

    source activate smrnaseq
    module load bowtie2/2.4.1

    [ ! -d ${name} ] && mkdir -p ${name}

    # bowtie
    bowtie2-build --quiet ${juncs} ${name}/${name} --threads ${task.cpus}
    """
}

process BOWTIE2_ALIGN_GENOME {

    label 'low'

    publishDir "$params.results/alignment/bam", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/alignment/bam", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results/alignment/genome_bed", mode : 'copy', pattern : "*.bed"
    publishDir "$params.results/alignment/genome_unmapped", mode : 'copy', pattern : "*.unmapped.fa"

    input : 
        val idx
        tuple val(sampleID), path(fasta)

    output : 
        tuple val(sampleID), path("*.unmapped.fa"), emit : unmapped_fa
        tuple val(sampleID), path("*.bed"), emit : genome_aligned_bed
        path("*.sorted.bam")
        path("*.sorted.bam.bai")
        path("*.bed")
        path("*.unmapped.fa")


    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    module load bowtie2/2.4.1

    name=\$(basename ${fasta} .fa)
    id=\$name.genome.aligned
    bam=\$id.sorted.bam
    bai=\$id.sorted.bam.bai
    bed=\$id.bed

    #############################################################################
    # Do Alignment with bowtie2 with base settings -a
    #############################################################################
    bowtie2 \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -a \\
        --un \$id.unmapped.fa \\
        -S > aligned.genome.sam 2> genome.log

    ########################################################################
    # Process alignment
    ########################################################################
    
    # sort sam file -> bam file
    samtools sort -m 1G -@ ${task.cpus} -o \$bam aligned.genome.sam

    # index bam file -> bam.bai 
    samtools index -@ ${task.cpus} \$bam
    
    # convert bam -> bed
    bam2bed < \$bam | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > \$bed
    """
}
    
process BOWTIE2_ALIGN_JUNCTION {

    label 'low'

    publishDir "$params.results/alignment/junc_bed", mode : 'copy', pattern : "*.bed"

    input : 
        val idx
        tuple val(sampleID), path(fasta)

    output : 
        tuple val(sampleID), path("*.bed"), emit : junction_aligned_bed

    script : 
    """
    #!/bin/bash

    source activate smrnaseq
    module load bowtie2/2.4.1

    name=\$(basename ${fasta} .fa)
    id=\$name.aligned.junctions
    bed=\$id.bed

    bowtie2 \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -a \\
        -S > aligned.junc.sam 2> junc.log

    # process alignment to junciton. bed -> junction bed (split the read across the junction)
    
    sam2bed < aligned.junc.sam  | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > tmp

    python3 ${params.bin}/junc_bed2bed.py -i tmp -o \$bed
    
    """
}