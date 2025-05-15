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

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
    'oras://community.wave.seqera.io/library/seqkit:2.10.0--9a5d37887d7c4e09' : 
    'community.wave.seqera.io/library/seqkit:2.10.0--03b4774218b4b7ef'}" 

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

// Classification with kraken
process KRAKEN2 {
    errorStrategy 'ignore'
    tag "classifying ${fastq_file.simpleName}"
    publishDir "${params.outdir}/kraken", mode:'copy'  

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/kraken2:2.14--508f454ac8eda8a9'
        : 'community.wave.seqera.io/library/kraken2:2.14--83aa57048e304f01'}"  

    input: 
        path kraken_db
        each path(fastq_file)  // [porechoped fastq]


    output:
        path "${filename}.classified.fastq"         , emit: fastq
        path "${filename}.pathogen_names.tsv"       , emit: tsv
        path "${filename}.report.txt"               , emit: txt

    script:
    filename = fastq_file.simpleName

    """
    kraken2 \\
	    --db $kraken_db  \\
	    --classified-out ${filename}.classified.fastq \\
	    --output ${filename}.pathogen_names.tsv  \\
	    --report ${filename}.report.txt  \\
	    --use-names   \\
	    --gzip-compressed $fastq_file
    """
}

// Get the tsv file from kraken
process LINUX_GREP {
    errorStrategy 'ignore'
    tag "generating final report"
    publishDir "${params.outdir}/grep", mode:'copy' 

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
    'oras://community.wave.seqera.io/library/grep:3.4--dc86f0f44bb51896' : 
    'community.wave.seqera.io/library/grep:3.4--845f7425f72cc8b4'}"

    input:
    path kraken_tsv

    output:
    path "${filename}.read_ids.txt"                 , emit: txt  // TSV file (may be empty)

    script:
    filename = kraken_tsv.simpleName

    """
    grep 'metapneumovirus' $kraken_tsv | cut -f 2 > ${filename}.read_ids.txt
    """
}

// Selected only hmpv virus reads
process SEQKIT_GREP {
    errorStrategy 'ignore'
    tag "generating final report"
    publishDir "${params.outdir}/seqkit_grep", mode:'copy' 

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
    'oras://community.wave.seqera.io/library/seqkit:2.10.0--9a5d37887d7c4e09' : 
    'community.wave.seqera.io/library/seqkit:2.10.0--03b4774218b4b7ef'}"
    
    input:
        path read_ids_txt
        path kraken_classified_fastq
        
    output:
        path "${filename}.hmpv_class.fastq.gz"       ,  emit: gz

    script:
    filename = kraken_classified_fastq.simpleName

    """
    seqkit grep \\
	    -f $read_ids_txt $kraken_classified_fastq | gzip  > ${filename}.hmpv_class.fastq.gz

    """
}

// Best reference selection
process AUTO_REF {
    errorStrategy 'ignore'
    tag "generating final report"
    publishDir "${params.outdir}/auto_ref", mode:'copy' 

    container "${workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ? 
    'docker://samordil/artic-multipurpose:1.6.2 : 
    'docker.io/samordil/artic-multipurpose:1.6.2'}"
    
    input:
        path multi_ref_fasta
        each path(fastq_file)
        

    output:
        path "${filename}.best_ref.fasta"       ,  emit: fasta
        path "*.json"                           ,  emit: json

    script:
    filename = fastq_file.simpleName

    """
   auto_ref.py \\
        --reads $fastq_file \\
	    --msa $multi_ref_fasta \\
	    --output ${filename}.best_ref.fasta \\
	    --threads $task.cpus \\
	    --min-coverage 80 \\
	    --preset map-ont
    """
}

process MINIMAP_SAMTOOLS {
    errorStrategy 'ignore'
    tag "generating final report"
    publishDir "${params.outdir}/consensus", mode:'copy' 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/minimap2_samtools:5919d63e7b60a09d' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

  tag {filename}

    input:
      path best_ref_fasta
      path fastq_gz   // .bam file

    output:
      path "${filename}.20X.consensus.fasta"      , emit: fasta

    script:
    filename = bam_file.simpleName

    """
    minimap2 \\
	-ax map-ont $best_ref_fasta $fastq_gz | samtools view -b -F 4 | samtools sort -o ${filename}.bam

    # Call the consensus
    samtools consensus --threads $task.cpus ${filename}.bam -aa -f fasta -d 20 -o ${filename}.20X.consensus.fasta    
    """
}
