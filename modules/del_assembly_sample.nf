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

// process HG_INDEX {
//     tag "Indexing human genome"
//     publishDir "${params.outdir}/index", mode:'copy'    

//     input:
//     // expecting -> [ path/dir/fasta.gz ]
//     path human_genome
    
//     output:
//     // Expected -> barcode01.fastq file
//     path "human_genome.mmi"     ,   emit: human_index_mmi 

//     script:

//     """
//     # Index
//     minimap2 -d human_genome.mmi $human_genome
//     """
// }

// STEP 00: Generate samplesheet
process GENERATE_SAMPLESHEET {
    tag "samplesheet generation"
   
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

// STEP 01: Generate index for the human genome
process INDEXING {
    tag "Indixing human genome"

    output:
        ind

        script:
        
        """
        minimap2 -d ref.mmi ref.fa                     # indexing
        """
}

// STEP 1: Remove adaptors using porechop
process REMOVE_HUMAN_READS {
    tag "Processing ${barcordeName}"
    publishDir "${params.outdir}/non_human_reads", mode:'copy'    

    input:
    // expecting -> [ path/dir/barcorde01 ]
    path human_genome_index
    each path(barcordDir)
    
    output:
    // Expected -> barcode01.fastq.gz file
    path "${barcordeName}.non_human_reads.fastq.gz"     ,   emit: fastq_gz 

    script:

    barcordeName = barcordDir.simpleName

    """
    # download human genome
    # wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
    
    # Concatenate the zipped fastq files
    cat ${barcordDir}/*.fastq.gz > ${barcordeName}.fastq.gz

    # Mapping
    minimap2 -ax map-ont $human_genome_index ${barcordeName}.fastq.gz > ${barcordeName}.sam

    # Filter out human reads
    samtools view -Sb ${barcordeName}.sam | samtools sort -o ${barcordeName}.sorted.bam
    samtools index ${barcordeName}.sorted.bam

    # Extract non-human reads
    samtools view -b -f 4 ${barcordeName}.sorted.bam > ${barcordeName}.non_human_reads.bam

    # Convert to fastq
    samtools fastq ${barcordeName}.non_human_reads.bam | gzip > ${barcordeName}.non_human_reads.fastq.gz

    """
}

// STEP 1: Remove adaptors using porechop
process PORECHOP {
    tag "Processing ${barcordeName}"
    publishDir "${params.outdir}/step1_porechop", mode:'copy'    

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

// STEP 2: Remove low quality reads using chopper
process CHOPPER { 
    tag "Processing ${barcordeName}"
    publishDir "${params.outdir}/step1_chopper", mode:'copy'    

     
    input:
    // expecting -> [ path/to/barcorde01.fastq.gz ] from porechop
    path porechop_fastq_gz
    
    output:
    // Expected -> barcode01.fastq.gz
    path "${barcordeName}.chopper.filtered.fastq.gz"    ,   emit: filtered_fastq_gz

    script:
    barcordeName = porechop_fastq_gz.simpleName

    """
    # adaptor removal 
        chopper \\
        --quality 9 \\
        --minlength 400 \\
        --threads $task.cpus \\
        --input $porechop_fastq_gz | gzip > ${barcordeName}.chopper.filtered.fastq.gz
    """
}

// STEP 3: Remove primer sequences 
process CUTADAPT { 
    tag "Processing ${barcordeName}"
    publishDir "${params.outdir}/step1_cutadapt", mode:'copy'    

    input:
    // expecting -> [ path/to/barcorde01 ]
    each path(chopped_fastq_gz)
    path foward_primers_txt
    path reverse_primers_txt
    
    output:
    // Expected -> barcode01.fastq file from chopper
    path "${barcordeName}.primer_free.fastq.gz"    ,   emit: trimmed_fastq_gz

    script:
    barcordeName = chopped_fastq_gz.simpleName

    """
    # primer removal 
    cutadapt \\
        --overlap 15 \\
        --cores 0 \\
        --minimum-length 50 \\
        -g file:$foward_primers_txt \\
        -a file:$reverse_primers_txt \\
        -o ${barcordeName}.primer_free.fastq.gz $chopped_fastq_gz
    """
}

// STEP 4: genome assembly 
// Assembly process
process FLYE { 
    tag "Processing ${barcordeName}"

    errorStrategy 'ignore'

    publishDir "${params.outdir}/step4_flye/", mode:'copy'    

     
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


// STEP 04: Remove adaptors using porechop
// process PORECHOP {
//     errorStrategy 'ignore'
//     tag "Processing ${filename}"
//     publishDir "${params.outdir}/step1_porechop", mode:'copy'    

//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/porechop:0.2.4--py39h2de1943_9' :
//         'biocontainers/porechop:0.2.4--py310h184ae93_9' }"

//     input:
//     // expecting -> [ path/dir/barcorde01 ]
//     path fastq_file
    
//     output:
//     // Expected -> barcode01.fastq file
//     path "${filename}.porechopped.fastq"     ,   emit: fastq

//     script:

//     filename = fastq_file.simpleName

//     """
//     # adaptor removal 
//      porechop \\
//         -i $fastq_file \\
//         --threads $task.cpus \\
//         --format fastq \\
//         -o ${filename}.porechopped.fastq
//     """
// }