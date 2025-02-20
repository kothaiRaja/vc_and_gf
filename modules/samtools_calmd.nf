process SAMTOOLS_CALMD {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"

    publishDir "${params.outdir}/calmd", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai),val(strandedness) 
          path(genome_fasta)
		  path (index)
		  
		  
    output:
    tuple val(sample_id), 
          path("${sample_id}_calmd.bam"), 
          path("${sample_id}_calmd.bam.bai"),
		  val(strandedness)
         

    script:
    """
    # Add NM and MD tags using samtools calmd
    samtools calmd -b ${bam} ${genome_fasta} > ${sample_id}_calmd.bam

    # Index the updated BAM file
    samtools index ${sample_id}_calmd.bam

    """
}