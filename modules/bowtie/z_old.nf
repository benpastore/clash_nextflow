process BOWTIE_INDEX {

    label 'low'
    
    //publishDir "$params.bowtie_index_path", mode : 'copy'
    publishDir "$params.index/bowtie", mode : 'copy'

    input :
        val genome
        val juncs
        val index
        val chrom_sizes
        val name

    output :
        val(index), emit : bowtie_index
        val(chrom_sizes), emit : bowtie_chrom_sizes
        path("${name}/*")

        //path("*.ebwt")
        //path("*chrom_sizes"), emit : bowtie_chrom_sizes
        //val("${params.index}/bowtie/${name}/${name}"), emit : bowtie_index

    script :
    """
    #!/bin/bash

    source activate smrnaseq

    [ ! -d ${name} ] && mkdir -p ${name}

    # bowtie
    bowtie-build --quiet ${genome} ${name}/${name} --threads ${task.cpus}
    bowtie-build --quiet ${juncs} ${name}/${name}_juncs --threads ${task.cpus}

    samtools faidx ${genome}
    cat ${genome}.fai | cut -f1,2 > ${name}/${name}_chrom_sizes
    rm ${genome}.fai
    """
}

process REMOVE_CONTAMINANT {
    
    label 'low'
    
    publishDir "$params.results/trim_galore/x_contaminant", mode : 'copy', pattern : "*.xc.fa"
    publishDir "$params.results/trim_galore/x_contaminant", mode : 'copy', pattern : "*.sam"

    input : 
        val idx 
        tuple val(sampleID), path(fasta)
    
    output : 
        tuple val(sampleID), path("*.xc.fa"), emit : xk_fasta
        path("*.sam")
        
    script :
    """
    #!/bin/bash
    
    source activate smrnaseq 
    
    name=\$(basename ${fasta} .fa)
    
    #build index 
    bowtie-build --quiet ${idx} reference --threads ${task.cpus}
    
    # align 
    bowtie \\
        -x reference \\
        -f ${fasta} \\
        --un \${name}.xc.fa \\
        -S > \${name}.aligned.contaminant.sam
    """

}

process BOWTIE_ALIGN_GENOME {

    label 'low'

    publishDir "$params.results/alignment/ntm", mode : 'copy', pattern : "*.ntm"
    publishDir "$params.results/alignment/rpm", mode : 'copy', pattern : "*.rpm"
    publishDir "$params.results/alignment/bw", mode : 'copy', pattern : "*.bw"
    publishDir "$params.results/alignment/bam", mode : 'copy', pattern : "*.bai"
    publishDir "$params.results/alignment/bam", mode : 'copy', pattern : "*.bam"
    publishDir "$params.results/alignment/rpkm", mode : 'copy', pattern : "*.rpkm"
    publishDir "$params.results/bowtie_unaligned", mode : 'copy', pattern : "*.unmapped.genome.junc.v0.m1.fa"
    publishDir "$params.results/logs", mode : 'copy', pattern : "*.log"

    input : 
        val idx
        tuple val(sampleID), path(fasta)
        val chrom_sizes

    output : 
        tuple val(sampleID), path("*.unmapped.genome.junc.v0.m1.fq"), emit : tailor_input
        tuple val(sampleID), path("*.unmapped.genome.junc.v0.m1.fa"), emit : unmapped_fa
        tuple val(sampleID), path("*.ntm"), emit : bowtie_alignment
        path("*.rpm"), emit : bowtie_rpm
        tuple val(sampleID), path("*.depth"), emit : normalization_constants
        path("*aligned*.log"), emit : alignment_logs
        path("*.sorted.bam")
        path("*.sorted.bam.bai")
        path("*.bw")
        path("*.rpkm")

    script : 
    mismatch_command = params.mismatch || params.mismatch == 0 ? "-v ${params.mismatch}" : ''
    multimap_command = params.multimap ? "-m ${params.multimap}" : ''
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    name=\$(basename ${fasta} .fa)
    id=\$name.aligned.v${mismatch}.m${multimap}
    bam=\$id.genome.sorted.bam
    bai=\$id.genome.sorted.bam.bai
    bed=\$id.bed
    ntm=\$id.ntm
    rpm=\$id.rpm
    rpkm=\$id.rpkm
    bw=\$id.bw
    log=\$id.log
    dep=\$id.depth

    #############################################################################
    # Do Alignment with -a --best --strata, get all aligned reads in best stratum
    #############################################################################
    bowtie \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -a \\
        --un \$name.tmp \\
        --best \\
        --strata \\
        ${mismatch_command} \\
        ${multimap_command} \\
        -S > aligned.genome.sam 2> genome.log
    
    bowtie \\
        -x ${idx}_juncs \\
        -f \$name.tmp \\
        -p ${task.cpus} \\
        -a \\
        --un \$name.unmapped.uni.fa \\
        --best \\
        --strata \\
        ${mismatch_command} \\
        ${multimap_command} \\
        -S > aligned.junc.sam 2> junc.log

    ##############################################################################
    # Do Alignment with -v 0, get reads that do not map perfectly, send to Tailor
    ##############################################################################
    bowtie \\
        -x ${idx} \\
        -f ${fasta} \\
        -p ${task.cpus} \\
        -v 0 \\
        -m 1 \\
        -a \\
        --un unmapped.genome.v0.tmp \\
        --best \\
        --strata \\
        -S > mapped.v0.sam
    python3 ${params.bin}/uniq_fasta_to_uniq_fastq.py unmapped.genome.v0.tmp > \$name.unmapped.genome.junc.v0.m1.fq

    #bowtie \\
    #    -x ${idx}_juncs \\
    #    -f unmapped.genome.v0.tmp \\
    #    -p ${task.cpus} \\
    #    -v 0 \\
    #    -m 1 \\
    #    -a \\
    #    --un \$name.unmapped.genome.junc.v0.m1.fa \\
    #    --best \\
    #    --strata \\
    #    -S > mapped.junc.v0.sam
    #python3 ${params.bin}/uniq_fasta_to_uniq_fastq.py \$name.unmapped.genome.junc.v0.m1.fa > \$name.unmapped.genome.junc.v0.m1.fq

    ########################################################################
    # Process alignment
    ########################################################################
    
    # sort sam file -> bam file
    samtools sort -m 1G -@ ${task.cpus} -o \$bam aligned.genome.sam

    # index bam file -> bam.bai 
    samtools index -@ ${task.cpus} \$bam
    
    # convert bam -> bed
    bam2bed < \$bam | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > aligned.genome.sorted.bed

    # process alignment to junciton. bed -> junction bed (split the read across the junction)
    sam2bed < aligned.junc.sam  | awk -F'\\t' -v OFS='\\t' '{split(\$4,a,":"); split(\$16,b,":"); print \$1,\$2,\$3,a[1],a[2],\$6,b[3] }' > tmp
    python3 ${params.bin}/junc_bed2bed.py -i tmp -o aligned.junc.bed

    # combine bed and junction bed -> bed
    cat aligned.genome.sorted.bed aligned.junc.bed > \$bed

    # bed -> ntm (split read count / number of locations mapped)
    python3 ${params.bin}/bed_to_ntm.py \$bed \$ntm

    # find depth from ntm file
    depth=\$(python3 ${params.bin}/add_columns.py \$ntm 4)
    echo \$depth > depth

    # normalize ntm file to rpm, remove multimappers and non-perfect alignments
    # awk '(\$7==1 && \$8==0)'
    echo -e "chrom\tstart\tend\tseq\tcount_rpm\tstrand" > header 
    cat \$ntm | awk -F'\\t' -v OFS='\\t' -v nTag=\$depth '{print \$1,\$2,\$3,\$4,1000000*(\$5/nTag),\$6}' > rpm.tmp
    cat header rpm.tmp > \$rpm

    # normalize ntm file to rpkm, remove multimappers and non-perfect alignments
    echo -e "chrom\tstart\tend\tseq\tcount_rpkm\tstrand" > header 
    cat \$rpm | awk -F'\\t' -v OFS='\\t' '{print \$1,\$2,\$3,\$4,\$5/length(\$4),\$6}' > rpkm.tmp
    cat header rpkm.tmp > \$rpkm

    # make bw file of all alignments
    cat rpm.tmp | awk -F'\\t' '{OFS="\\t"; if (\$2>=\$3) print \$1,\$2,\$3+1,\$4,\$5; else print \$1,\$2,\$3,\$4,\$5 }' | sort-bed - > tmp.sorted
    bedops --partition tmp.sorted | bedmap --echo --sum --delim '\\t' - tmp.sorted | awk -F '\\t' '{OFS="\\t"; print \$1,\$2,\$3,1*\$4}' > tmp.bg
    ${params.bin}/bedGraphToBigWig tmp.bg ${chrom_sizes} \$bw

    # find total reads sequenced
    total_reads=\$(cat ${fasta} | grep "^>" | sed -e 's/.*://g' | ${params.bin}/addCols stdin)
    
    # compute total aligned reads (decimal percent)
    percent_aligned=\$(python3 -c "print(\$depth/\$total_reads)")
    
    # Alignment summary
    echo -e "sample\ttotal_reads\tdepth\tpercent_aligned" > \$log
    echo -e \$id'\t'\$total_reads'\t'\$depth'\t'\$percent_aligned >> \$log

    # Normalization constants 
    echo -e 'total\t'\$depth > \$dep

    # fin
    """
}

// Make Counts Master Table
process BOWTIE_ALIGNMENT_MASTER_TABLE {

    label 'medium'

    publishDir "$params.results/alignment/alignment_master_table", mode : 'copy', pattern : "*tsv"

    input :
        val project_name
        val counts

    output : 
        path("*count*tsv"), emit : tables
        path("*tsv"), emit : unnormalized_master_table_ch
    
    script :
    mismatch = params.mismatch || params.mismatch == 0 ? "${params.mismatch}" : ''
    multimap = params.multimap ? "${params.multimap}" : ''
    """
    #!/bin/bash

    source activate smrnaseq

    counts=${project_name}.bowtie.aligned.v${mismatch}.m${multimap}

    time python3 ${params.bin}/make_master.py -f "${counts}" -o \$counts
    """
}