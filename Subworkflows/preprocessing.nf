nextflow.enable.dsl = 2

// Include required processes
include { CONCAT_FASTQ } from '../modules/concat_fastq.nf'
include { FASTQC_RAW } from '../modules/fastqc_raw.nf'
include { TRIM_READS } from '../modules/trim_reads.nf'
include { MultiQC_quality } from '../modules/multiqc_quality.nf'

workflow PREPROCESSING {
    take:
    samplesheet

    main:

    log.info "🛠 Starting Preprocessing Steps..."

    // Step 1: Read samplesheet
    samples_ch = Channel.fromPath(samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)], row.strandedness) }

    // Step 2: Concatenate FASTQ files if required
    concatenated_reads_ch = params.concatenate ? 
        CONCAT_FASTQ(samples_ch) : 
        samples_ch.map { sample_id, reads, strandedness -> tuple(sample_id, reads[0], reads[1], strandedness) }

    // Step 3: Perform FastQC on raw reads
    qc_results_ch = FASTQC_RAW(concatenated_reads_ch)

    // Step 4: Trim reads using Fastp
    trimmed_reads_ch = TRIM_READS(concatenated_reads_ch)

    // Step 5: Collect QC reports for MultiQC
    qc_files_ch = qc_results_ch.map { it[1] + it[2] }.flatten()
    fastp_files_ch = trimmed_reads_ch.fastp_reports.map { it[1..-1] }.flatten()
    combined_channel = qc_files_ch.concat(fastp_files_ch).collect()
    multiqc_quality = MultiQC_quality(combined_channel)

    log.info "✅ Preprocessing Completed."

    emit:
        trimmed_reads    = trimmed_reads_ch.trimmed_reads  // Ensuring a single-channel output
        fastp_reports    = trimmed_reads_ch.fastp_reports  // Keeping fastp reports separate
        qc_reports       = qc_results_ch                   // Ensuring QC reports are properly outputted
		multiqc        = multiqc_quality
}
