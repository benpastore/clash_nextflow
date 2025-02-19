#!/usr/bin/env nextflow

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

/*
////////////////////////////////////////////////////////////////////
set path to bin and index
////////////////////////////////////////////////////////////////////
*/
params.bin = "${params.base}/bin"
params.index = "${params.base}/index"
if (params.outdir)   { ; } else { exit 1, 'Results path not specified!' }
params.date = new Date().format( 'yyyyMMdd' )
params.results = "${params.outdir}/${params.date}"
params.outprefix = "${params.date}"

/*
////////////////////////////////////////////////////////////////////
Included modules
////////////////////////////////////////////////////////////////////
*/

include { TRIM_GALORE } from './modules/trimgalore/main.nf'
include { TAILOR_INDEX } from './modules/tailor/main.nf'
include { TAILOR_ALIGN } from './modules/tailor/main.nf'
include { BWA_INDEX } from './modules/bwa/main.nf'
include { BWA_MEM } from './modules/bwa/main.nf'
include { PROCESS_BWA_MEM } from './modules/bwa/main.nf'
include { BOWTIE_INDEX_GENOME } from './modules/bowtie/main.nf'
include { BOWTIE_ALIGN_GENOME } from './modules/bowtie/main.nf'
include { MERGE_TAILOR_BWA } from './modules/merge_alignments/main.nf'
include { RNADUPLEX } from './modules/RNAduplex/main.nf'

/*
////////////////////////////////////////////////////////////////////
Workflows
////////////////////////////////////////////////////////////////////
*/

workflow trim_galore {

    take : 
        data 

    main : 
        TRIM_GALORE( data )
    
    emit : 
        trim_galore_output = TRIM_GALORE.out.trimmed_fq

}

workflow bowtie_align {

    take : 
        genome
        reads 
        index

    main : 
        genome_fasta = file("${genome}")
        genome_name = "${genome_fasta.baseName}"
        bowtie_index = "${index}/bowtie/${genome_name}"
        bowtie_exists = file(bowtie_index).exists()

        if ( bowtie_exists ){
                index_ch = "${bowtie_index}/${genome_name}"
                chrom_sizes = "${bowtie_index}/${genome_name}_chrom_sizes"

        } else {
            BOWTIE_INDEX_GENOME( 
                genome, 
                "${bowtie_index}/${genome_name}", 
                "${bowtie_index}/${genome_name}_chrom_sizes", 
                genome_name)

            index_ch = BOWTIE_INDEX_GENOME.out.bowtie_index
            chrom_sizes = BOWTIE_INDEX_GENOME.out.bowtie_chrom_sizes
        }

        BOWTIE_ALIGN_GENOME(index_ch, reads)

    emit : 
        unaligned = BOWTIE_ALIGN_GENOME.out.unmapped

}

workflow bwa_align {

    take :
        genome
        reads 
        index

    main :
        genome_fasta = file("${genome}")
        genome_name = "${genome_fasta.baseName}"
        bwa_index_path = "${index}/bwa/${genome_name}"
        bwa_index = "${bwa_index_path}/${genome_name}"
        bwa_chrom_sizes = "${bwa_index_path}/chrom_sizes.txt"

        bwa_exists = file(bwa_index_path).exists()

        if (bwa_exists == true){
            bwa_build = false
        } else {
            bwa_build = true
        }

        if (bwa_build == true){
            BWA_INDEX(
                genome,
                "${bwa_index}", 
                genome_name )

            index_ch = BWA_INDEX.out.bwa_index
        } else {
            index_ch = "${bwa_index}"
        }

        BWA_MEM( index_ch, reads )

    emit :
        sam = BWA_MEM.out.bwa_sam
}

workflow tailor_align {

    take : 
        genome
        reads
        index 

    main : 
        genome_fasta = file("${genome}")
        genome_name = "${genome_fasta.baseName}"
        index_path = "${index}/tailor/${genome_name}"
        index = "${index_path}/${genome_name}"
        chrom_sizes = "${index_path}/chrom_sizes.txt"

        exists = file(index_path).exists()

        if (exists == true){
            build = false
        } else {
            build = true
        }

        if (build == true){
            TAILOR_INDEX( genome, index, genome_name )
            index_ch = TAILOR_INDEX.out.tailor_index
        } else {
            index_ch = index 
        }

        TAILOR_ALIGN( index_ch, reads)
    
    emit : 
        bed = TAILOR_ALIGN.out.tailor_bed

}

workflow merge_tailor_bwa {
    
    take : 
        bwa
        tailor
    
    main : 
        MERGE_TAILOR_TO_BWA( bwa, tailor )
    
    emit : 
        merge_tsv = MERGE_TAILOR_BWA.out.merge_tsv


}

workflow rnaduples { 

    take : 
        data
        transcripts
        smrnas
    
    main : 
        RNADUPLEX( data, transcripts, smrnas )

    emit : 
        tsv = RNADUPLEX.out.tsv 

}


/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    /*
    ////////////////////////////////////////////////////////////////////
    Validate mandatory inputs (design, genome, junctions, results, outprefix)
    ////////////////////////////////////////////////////////////////////
    */
    if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
    if (params.targets)   { ch_targets = file(params.targets, checkIfExists: true) } else { exit 1, 'Target fasta not specified!' }
    if (params.smrnas)   { ch_smrnas = file(params.smrnas, checkIfExists: true) } else { exit 1, 'Target fasta not specified!' }
    if (params.outprefix) { ; } else {'Outprefix not specified! Defaulting to smRNA_analysis'; params.outprefix = 'smRNA_analysis' }


    sample_ch = Channel.fromPath(params.reads)
                .map { file -> 
                    def sampleID = file.getSimpleName()
                    return [sampleID, file]
                }

    sample_ch.view()
    

    // trimming and qc 
    if ( params.format == 'fastq' & params.trimming ) {
        trim_galore( sample_ch )
        fastas = trim_galore.out.trim_galore_output
    } else {
        fastas = sample_ch
    }

    // Step 1. align to genome and keep reads that do not align
    bowtie_align( params.genome, fastas, params.index)

    unal = bowtie_align.out.unaligned

    // Step 2. align to targets with bwa, process with pysam script 
    bwa_align(params.targets, unal, params.index)

    PROCESS_BWA_MEM( bwa_align.out.sam )

    // Step 3. take portion that did not align and align it to smRNA reference using tailor
    tailor_align(params.smrnas, PROCESS_BWA_MEM.out.bwa_bed, params.index)

    // Step 4. merge tailor alignments to bwa alignment
    merge_tailor_bwa( PROCESS_BWA_MEM.out.bwa_bed, tailor_align.out.bed )

    // Step 5. RNAduplex
    rnaduplex( merge_tailor_bwa.out.merge_tsv, params.targets, params.smrnas )

    


}

