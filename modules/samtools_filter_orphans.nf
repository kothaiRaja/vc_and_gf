process SAMTOOLS_FILTER_ORPHANS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/filtered_bam", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai),val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai"),val(strandedness)

    script:
    """
    # Filter out orphan reads (retain properly paired reads)
    samtools view -h -f 2 ${sorted_bam} | samtools view -b - > ${sample_id}_filtered.bam

    # Index the filtered BAM file
    samtools index ${sample_id}_filtered.bam
    """
}
