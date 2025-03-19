#!/usr/bin/env nextflow

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

/*
========================================================================================
    A Nextflow process block. Process names are written, by convention, in uppercase.
    This convention is used to enhance workflow readability.

    Tutorials
    https://timkahlke.github.io/LongRead_tutorials/ 
    https://colauttilab.github.io/NGS/deNovoTutorial.html
    https://star-protocols.cell.com/protocols/1799
========================================================================================
*/

// STEP 01: Generate samplesheet
process GENERATE_SAMPLESHEET {
    errorStrategy 'ignore'
    tag "samplesheet generation"

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
    'docker://samordil/python-pandas-dateutil:1.0.0' : 
    'docker.io/samordil/python-pandas-dateutil:1.0.0'}"
    
    input:
        path fastq_dir
        path metadata_tsv

    output:
        path "samplesheet.csv"                          , emit: samplesheet
        path "*without*.csv"    , optional: true        , emit: csv


    script:     // This script is bundled with the pipeline, in kwtrp-peo/viralphyl/bin/
    def metadata  = metadata_tsv ? "--metadata $metadata_tsv" : ""

    """
   generate_samplesheet.py \\
        --directory $fastq_dir \\
        $metadata \\
        --output samplesheet.csv
    """
}

// STEP 02: Indexing the human genome
process HG_INDEXING {
    tag "Indexing human genome reference"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }" 

    input: 
        path ref_genome

    output:
        path "indexed_human_ref.mmi"           , emit: indexed_file
    script:
    def human_genome  = ref_genome ? "$ref_genome" : ""

    """
    
    if [[ -n "$human_genome" ]]; then
        # Index the human genome
        minimap2 -d indexed_human_ref.mmi $human_genome
    else
        # download the genome
        wget -O human.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
        
        # Index the genome
        minimap2 -d indexed_human_ref.mmi human.gz
    fi

    """
}

// STEP 03: remove human reads
process REMOVE_HUMAN_READS {
    // errorStrategy 'ignore'
    tag "Processing ${sample_id}"
    publishDir "${params.outdir}/non_human_reads", mode:'copy'  

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"  

    input:
    path indexed_ref // [ index ]
    tuple(val(sample_id), path(fastq_dir)) // expecting [ [sample_id], [fastq_dir] ]
    
    output:
    // Expected -> barcode01.fastq.gz file
    path "${sample_id}.non_human_reads.fastq.gz"     ,   emit: fastq_gz 

    script:
    
    """
    # Concatenate the zipped fastq files
    cat ${fastq_dir}/*.fastq.gz > ${sample_id}.fastq.gz

    minimap2 -ax map-ont $indexed_ref ${sample_id}.fastq.gz | \
        samtools view -b -f 4 | \
        samtools sort | samtools fastq | gzip > ${sample_id}.non_human_reads.fastq.gz
    """
}

// STEP 04: Remove adaptors using porechop
process PORECHOP {
    tag "Processing ${filename}"
    publishDir "${params.outdir}/step1_porechop", mode:'copy'    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h2de1943_9' :
        'biocontainers/porechop:0.2.4--py310h184ae93_9' }"

    input:
    // expecting -> [ path/dir/barcorde01 ]
    path fastq_file
    
    output:
    // Expected -> barcode01.fastq file
    path "${filename}.porechopped.fastq"     ,   emit: fastq

    script:

    filename = fastq_file.simpleName

    """
    # adaptor removal 
     porechop \\
        -i $fastq_file \\
        --threads $task.cpus \\
        --format fastq \\
        -o ${filename}.porechopped.fastq
    """
}

// STEP 04: taxonomic classification  : Downloading mash database in not provided
process MASH_CLASSIFICATION {
    errorStrategy 'ignore'
    tag "classifying ${fastq_file.simpleName}"
    publishDir "${params.outdir}/classify", mode:'copy'  

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--h7c2a333_8' :
        'biocontainers/mash:2.3--hb105d93_9' }"

    input: 
        path mash_database
        each(path fastq_file)   // [fastq_dir]


    output:
        path "${filename}.mash.results.txt"         , emit: txt

    script:
    def mash_db  = mash_database ? "$mash_database" : ""
    filename = fastq_file.simpleName

    """
    
    if [[ -n "$mash_db" ]]; then
        # Create a Mash Sketch from FASTQ File
        # mash sketch -o ${filename}.msh $fastq_file

        # Run mash screen Using the Prebuilt Database
        mash screen -w -p 16 $mash_db $fastq_file | sort -gr > ${filename}.mash.results.txt
    else
        # download the genome
        wget -O mash_downloaded_db.msh https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21.s1000.msh

        # # create sketch from the sample
        # mash sketch -o ${filename}.msh $fastq_file

        # Run mash screen Using the Prebuilt Database
        mash screen -w -p 16 mash_downloaded_db $fastq_file | sort -gr > ${filename}.mash.results.txt
    fi

    """
}
