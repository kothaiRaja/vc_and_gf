process TRIM_READS {
    tag { sample_id }
    container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--h125f33a_5"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_fastp.*"

    input:
    tuple val(sample_id), path(r1), path(r2), val(strandedness)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), val(strandedness), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports

    script:
    """
    fastp -i ${r1} -I ${r2} \
      -o trimmed_${sample_id}_R1.fastq.gz \
      -O trimmed_${sample_id}_R2.fastq.gz \
      --detect_adapter_for_pe \
      --adapter_sequence auto \
      --adapter_sequence_r2 auto \
      --length_required 50 \
      --cut_front --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --html ${sample_id}_fastp.html \
      --json ${sample_id}_fastp.json
    """
}