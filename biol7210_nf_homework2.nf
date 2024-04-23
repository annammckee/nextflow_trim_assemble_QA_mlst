#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Set parameter values
params.reads = "$baseDir/data/*_{R1,R2}.fastq.gz"
params.outdir = "$baseDir/results"

// Process that trims the input FastQ files using Trimmomatic
process TrimReads {
    tag "${pair_id}" // Tagging helps in identifying the job by the pair_id in the log files
    publishDir "${params.outdir}/trimmed", mode: 'copy'// Specifies output directory for trimmed files

    input:
    tuple val(pair_id), path(read1), path(read2) // Inputs are tuples of pair_id and paths for R1 and R2

    output:
    tuple val(pair_id), path("${pair_id}_R1.trimmed.fastq.gz"), path("${pair_id}_R2.trimmed.fastq.gz") // Outputs trimmed FastQ files

    script:
    """
    trimmomatic PE -phred33 \\
        $read1 $read2 \\
        ${pair_id}_R1.trimmed.fastq.gz ${pair_id}_R1.unpaired.fastq.gz \\
        ${pair_id}_R2.trimmed.fastq.gz ${pair_id}_R2.unpaired.fastq.gz \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:36
    """// trimmomatic command for quality trimming and adapter removal
}

// Process for assembling the trimmed reads into contigs using SPAdes
process AssembleFasta {
    tag "${pair_id}" // Tags jobs with the pair_id
    publishDir "${params.outdir}/assembly", mode: 'copy' // Specifies output directory for assemblies

    input:
    tuple val(pair_id), path(r1_trimmed), path(r2_trimmed) // Inputs are the trimmed read files

    output:
    tuple val(pair_id), path("${pair_id}.fasta"), emit: fasta_file // Outputs assembled fasta files and emit them to a channel

    script:
    """
    spades.py -1 $r1_trimmed -2 $r2_trimmed -o ${pair_id}.assembly
    mv ${pair_id}.assembly/contigs.fasta ${pair_id}.fasta
    """ // SPAdes assembly command followed by moving the result to a named fasta file
}

// Quality assessment using QUAST
process QA {
    tag "${pair_id}" // Tags jobs with the pair_id
    publishDir "${params.outdir}/qa", mode: 'copy' // Specifies output directory for QUAST reports

    input:
    tuple val(pair_id), path(fasta_file) // Input fasta files for quality assessment

    output:
    path ("${pair_id}_quast") // Output directory for QUAST reports

    script:
    """
    quast ${fasta_file} -o ${pair_id}_quast
    """  // QUAST command to evaluate assembly quality
}

// Multi-Locus Sequence Typing (MLST)
process MLST {
    tag "${pair_id}" // Tags jobs with the pair_id
    publishDir "${params.outdir}/mlst", mode: 'copy' // Specifies output directory for MLST results

    input:
    tuple val(pair_id), path(fasta_file)  // Input fasta file for MLST analysis

    output:
    path "${pair_id}_mlst.csv" // Output MLST results in CSV format

    script:
    """
    mlst --csv ${fasta_file} > ${pair_id}_mlst.csv
    """ // MLST command that generates CSV output
}

// Defines the workflow block that orchestrates the execution of the defined processes
workflow {
    read_pairs = Channel.fromFilePairs(params.reads, size: 2, flat: true) // Channel to pair R1 and R2 reads

    trimmed_reads = read_pairs
        | TrimReads // Trimming read pairs

    assembled_fasta = trimmed_reads
        | AssembleFasta // Assembling trimmed reads into fasta files

    qa_results = assembled_fasta
        | QA // Running quality assessment on assembled fasta files

    mlst_results = assembled_fasta
        | MLST // Running MLST on assembled fasta files

    // Output a completion message for each assembled fasta file
    assembled_fasta
        .view { it -> "Assembly completed: ${it}" }
}
