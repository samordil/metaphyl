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
            GENERATE_FINAL_REPORT;
            PORECHOP;
            FLYE;
            QUAST;
            MULTIQC;
            SEQKIT;
            KRAKEN2;
            AUTO_REF;
            LINUX_GREP;
            SEQKIT_GREP;
            MINIMAP_SAMTOOLS
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

    // Define kraken directory channel
    Channel                                                     // Get kraken db directory
        .fromPath(params.kraken_db, type: 'dir')
        .set { ch_kraken_db }

    Channel                                                     // Get multiref file
        .fromPath(params.multi_ref_fasta)
        .set { ch_multiref }
 
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
    PORECHOP (REMOVE_HUMAN_READS.out.fastq)

    // MODULE: Denove genome assembly
    FLYE (PORECHOP.out.fastq_gz)

    // Generate assembly statistics
    QUAST(
      FLYE.out.assembly_fasta
    )

    // Aggregage report for visualization
    MULTIQC (
      QUAST.out.quast_report_dir.collect(),
      params.multiqc_title
    )

    // Get the longest contigs per sample
    SEQKIT (
      FLYE.out.assembly_fasta
    )

    // MODULE: Classify reads
    MASH_CLASSIFICATION ( ch_mash_db,
                        REMOVE_HUMAN_READS.out.fastq
                        )

    // MODULE: Aggregate mash txts
    GENERATE_FINAL_REPORT (
        MASH_CLASSIFICATION.out.txt.collect()
    )

    /*
    Classification using kraken2
    */
    // 
    // Classification using kraken
    KRAKEN2 (
        ch_kraken_db,
        PORECHOP.out.fastq_gz
    )

    // Get the read ids of hmpv
    LINUX_GREP (
        KRAKEN2.out.tsv
    )

    // Only pass TSVs with ≥1000 line to SEQKIT
    LINUX_GREP.out.txt
        .filter { file -> file.countLines() > 1000 }  // Skip empty files
        .map { tuple(it.getSimpleName(), it) }
        .set { ch_read_ids }

    KRAKEN2.out.fastq
        .map { tuple(it.getSimpleName(), it) }
        .set { ch_kraken2_fastq}

    // Join operation that will only keep matched pairs
    ch_read_ids.join(ch_kraken2_fastq).set { ch_paired_samples }

        
    // Extract hmpv reads only
    SEQKIT_GREP (
        ch_paired_samples
    )

    // SEQKIT_GREP.out.gz.view()

    // Get the best reference
    AUTO_REF (
        ch_multiref,
        SEQKIT_GREP.out.gz
    )

    // Prepare channel for MINIMAP_SAMTOOLS process
    AUTO_REF.out.fasta
        .map { tuple(it.getSimpleName(), it) }
        .set { ch_best_ref }

    SEQKIT_GREP.out.gz
        .map { tuple(it.getSimpleName(), it) }
        .set { ch_seqkit_grep_gz }

    // Combine best reference with the corresponding fastq
    ch_best_ref.join(ch_seqkit_grep_gz).set { ch_best_ref_fastq }

    // Generate consensus
    MINIMAP_SAMTOOLS (
        ch_best_ref_fastq
    )
}

