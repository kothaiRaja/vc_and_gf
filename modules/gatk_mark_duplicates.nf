process GATK_MARK_DUPLICATES {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_dup_metrics.txt"

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bam_index),val(strandedness)

    output:
    tuple val(sample_id), 
          path("${sample_id}_marked_duplicates.bam"), 
          path("${sample_id}_marked_duplicates.bai"),
		  val(strandedness),
          path("${sample_id}_dup_metrics.txt")

    script:
    """
    gatk MarkDuplicates \
        -I ${sorted_bam} \
        -O ${sample_id}_marked_duplicates.bam \
        -M ${sample_id}_dup_metrics.txt \
        --CREATE_INDEX true \
		--REMOVE_DUPLICATES ${params.remove_duplicates ? 'true' : 'false'} \
        --VALIDATION_STRINGENCY ${params.validation_stringency ?: 'LENIENT'}
		
	# Check if output BAM is generated
    if [ ! -s ${sample_id}_marked_duplicates.bam ]; then
        echo "Error: Marked duplicates BAM file not generated for ${sample_id}" >&2
        exit 1
    fi

    # Check if metrics file is generated
    if [ ! -s ${sample_id}_dup_metrics.txt ]; then
        echo "Error: Duplicate metrics file not generated for ${sample_id}" >&2
        exit 1
    fi

    """
}