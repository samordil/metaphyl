#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
 include{   GENERATE_SAMPLESHEET;
            HG_INDEXING;
            REMOVE_HUMAN_READS;
            MASH_CLASSIFICATION;
            GENERATE_FINAL_REPORT
          } from "./modules/assembly_processes"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
    nextflow run main.nf --fastq_dir data \
        --outdir Results \
        --human_genome data/GCF_000001405.39_GRCh38.p13_genomic.fna.gz \
        -profile local --mash_db data/mash_local_db.msh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {

    // Define input channel for an optional tsv metadata file
    if (params.metadata_tsv) { ch_metadata_tsv = file(params.metadata_tsv) } else { ch_metadata_tsv = [] }

    // Define human reference channel for an optional reference genome
    if (params.human_genome) { ch_human_genome = file(params.human_genome) } else { ch_human_genome = [] }

    // Define human reference channel for an optional reference genome
    if (params.mash_db) { ch_mash_db= file(params.mash_db) } else { ch_mash_db = [] }
 
    // Define fastq_pass directory channel
    Channel                                                     // Get raw fastq directory
        .fromPath(params.fastq_dir, type: 'dir', maxDepth: 1)
        .set { ch_fastq_data_dir }

    // MODULE: Run bin/viraphly_samplesheet_generator.py to generate samplesheet
    GENERATE_SAMPLESHEET (
        ch_fastq_data_dir,               // raw reads directory channel
        ch_metadata_tsv                  // tsv metadata channel (can be an empty channel)
    )

     GENERATE_SAMPLESHEET.out.samplesheet
        .splitCsv(header:true)
        .map { row -> tuple(row.strain_id, file(row.fastq_dir)) }
        .set { ch_samplesheet }

    // MODULE: Index the human reference genome
    HG_INDEXING ( ch_human_genome )

    // MODULE: Remove human reads
    REMOVE_HUMAN_READS(
        HG_INDEXING.out.indexed_file.map {file(it)},
        ch_samplesheet
    )

    // MODULE: Remove adapters
    // PORECHOP (REMOVE_HUMAN_READS.out.fastq_gz)

    // MODULE: Classify reads
    MASH_CLASSIFICATION ( ch_mash_db,
                        REMOVE_HUMAN_READS.out.fastq
                        )

    // MODULE: Aggregate mash txts
    GENERATE_FINAL_REPORT (
        MASH_CLASSIFICATION.out.txt.collect()
    )
}
