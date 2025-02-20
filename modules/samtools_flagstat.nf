process SAMTOOLS_FLAGSTAT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_flagstat.*"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai), val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt"), path("${sample_id}_stats_report.txt"), val(strandedness)

    script:
    """
    # Validate BAM and index file
    if [ ! -s ${sorted_bam} ]; then
        echo "Error: BAM file is empty or missing for ${sample_id}" >&2
        exit 1
    fi
    
    if [ ! -s ${bai} ]; then
        echo "Error: BAM index (.bai) file is missing for ${sample_id}" >&2
        exit 1
    fi

    # Validate BAM file to ensure it is not corrupted
    samtools quickcheck -v ${sorted_bam} || (echo "BAM file validation failed for ${sample_id}" && exit 1)

    # Run samtools flagstat to generate alignment metrics
    samtools flagstat ${sorted_bam} > ${sample_id}_flagstat.txt

    # Generate additional quality metrics with samtools stats
    samtools stats ${sorted_bam} > ${sample_id}_stats_report.txt
    """
}
