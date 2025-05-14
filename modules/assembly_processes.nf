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
        wget -O GCF_000001405.39_GRCh38.p13_genomic.fna.gzftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
        
        # Index the genome
        minimap2 -d indexed_human_ref.mmi GCF_000001405.39_GRCh38.p13_genomic.fna.gz
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
    path "${sample_id}.non_human_reads.fastq"     ,   emit: fastq

    script:
    
    """
    # Concatenate the zipped fastq files
    cat ${fastq_dir}/*.fastq.gz > ${sample_id}.fastq.gz

    minimap2 -ax map-ont $indexed_ref ${sample_id}.fastq.gz | \
        samtools view -b -f 4 | \
        samtools sort | samtools fastq > ${sample_id}.non_human_reads.fastq
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
        each path(fastq_file)  // [fastq_dir]


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
        wget -O refseq.genomes.k21s1000.msh https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh

        # # create sketch from the sample
        # mash sketch -o ${filename}.msh $fastq_file

        # Run mash screen Using the Prebuilt Database
        mash screen -w -p 16 refseq.genomes.k21s1000.msh $fastq_file | sort -gr > ${filename}.mash.results.txt
    fi

    """
}


// STEP 05: Generate the concatenated tsv file
process GENERATE_FINAL_REPORT {
    errorStrategy 'ignore'
    tag "generating final report"
    publishDir "${params.outdir}/final_mash", mode:'copy' 

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
    'docker://samordil/python-pandas-dateutil:1.0.0' : 
    'docker.io/samordil/python-pandas-dateutil:1.0.0'}"
    
    input:
        path tsv_file

    output:
        path "mash_summary_out.tsv"   ,      emit: tsv

    script:     // This script is bundled with the pipeline, in kwtrp-peo/viralphyl/bin/

    """
   mash_summary.py \\
        --input $tsv_file \\
        --output mash_summary_out.tsv
    """
}


// Remove adaptors using porechop
process PORECHOP {
    tag "Processing ${barcordeName}"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/porechop", mode:'copy'    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h2de1943_9' :
        'biocontainers/porechop:0.2.4--py310h184ae93_9' }"

    input:
    // expecting -> [ path/dir/barcorde01 ]
    path barcordDir
    
    output:
    // Expected -> barcode01.fastq file
    path "${barcordeName}.porechopped.fastq.gz"     ,   emit: fastq_gz 

    script:

    barcordeName = barcordDir.simpleName

    """
    # adaptor removal 
     porechop \\
        -i $barcordDir \\
        --threads $task.cpus \\
        --format fastq.gz \\
        -o ${barcordeName}.porechopped.fastq.gz
    """
}

// genome assembly 
process FLYE { 
    tag "Processing ${barcordeName}"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/flye/", mode:'copy'    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye%3A2.9.6--py39h475c85d_0' :
        'biocontainers/flye:2.9.6--py39h475c85d_0' }"

    input:
    // expecting -> [ path/to/barcorde01.fastq.gz ] from cutadapt
    path trimmed_fasta_gz
    
    output:
    // Expected -> RSVB_barcode01.fastq file
    path "${barcordeName}.assembly.fasta"       , emit: assembly_fasta

    script:

    barcordeName = trimmed_fasta_gz.simpleName

    """
    flye \\
        --nano-raw  $trimmed_fasta_gz \\
        --genome-size 13k \\
        --threads $task.cpus \\
        --iterations 3 \\
        --out-dir $barcordeName

    # Rename the assembly.fasta file to the sample barcode name
    mv ${barcordeName}/assembly.fasta ${barcordeName}.assembly.fasta
    """
}

// Get assembly statistics using quast
process QUAST { 
    tag "Processing ${barcordeName}"
    publishDir "${params.outdir}/quast/", mode:'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a5/a515d04307ea3e0178af75132105cd36c87d0116c6f9daecf81650b973e870fd/data' :
        'community.wave.seqera.io/library/quast:5.3.0--755a216045b6dbdd' }"

    input:
    // expecting -> [ path/to/assembly_fasta ] from flye
    path assembly_fasta
    
    output:
    // Expected -> quast_result dir for multiqc
    path "${barcordeName}"      , emit: quast_report_dir

    script:
    barcordeName = assembly_fasta.simpleName

    """
    quast.py \\
    --threads $task.cpus \\
    --no-plots \\
    --no-html \\
    --no-icarus \\
    --no-snps \\
    --no-read-stats \\
    --output-dir $barcordeName $assembly_fasta 

    # mv ${barcordeName}/report.tsv ${barcordeName}.report.tsv
    """
}

// Get assembly statistics using quast
process MULTIQC { 
    tag "Aggregating all reports"
    publishDir "${params.outdir}/multiqc/", mode:'copy'    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.28--pyhdfd78af_0' :
        'biocontainers/multiqc:1.28--pyhdfd78af_0' }"

    input:
    // expecting -> [ path/to/assembly_fasta ] from flye
    path all_reports
    val title_name
    
    output:
    // Expected -> quast_result dir for multiqc
    path "flye_assembly_report.html"      , emit: multqc_report

    script:
    """
    multiqc $all_reports \\
    --title "$title_name" \\
    --filename flye_assembly_report.html

    """
}

// Get the longest contigs
process SEQKIT { 
    tag "Processing ${barcordeName}"
    publishDir "${params.outdir}/longest_contigs", mode:'copy'  

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"  

    input:
    // expecting -> [ path/to/barcorde01 ]
    path assembly_fasta
    
    output:
    // Expected -> RSVB_barcode01.fastq file from porechop
    path "${barcordeName}_longest_contig.fasta"       , emit: fasta

    script:
    barcordeName = assembly_fasta.simpleName

    """
    # Get the longest contig
    seqkit sort -l \\
        -r $assembly_fasta | seqkit head -n 1 > ${barcordeName}_longest_contig.fasta
    """
}