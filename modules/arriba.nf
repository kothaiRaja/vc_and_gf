process ARRIBA {
    tag { sample_id }
    
    container "https://depot.galaxyproject.org/singularity/arriba%3A2.4.0--hdbdd923_3"
    publishDir "${params.outdir}/ARRIBA", mode: 'copy'

    input:
	tuple val(sample_id),path(log_final), val(strandedness)
	tuple val(sample_id),path(log_out), val(strandedness)
	tuple val(sample_id),path(bam), val(strandedness)
	tuple val(sample_id),path(chimeric_sam), val(strandedness)
	tuple val(sample_id), path(log_progress) , val(strandedness)
    tuple val(sample_id), path(splice_junctions), val(strandedness) 
    path fasta                                     
    path gtf                                       
    path blacklist                                 
    path known_fusions  
                           

    output:
	tuple val(sample_id), path("${sample_id}.fusions.tsv"), val(strandedness), emit: fusions
	tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), val(strandedness), emit: fusions_discarded


    script:
    """
    arriba \
         -x $bam \
         -c ${sample_id}_Chimeric.out.sam \
         -a $fasta \
         -g $gtf \
         -b $blacklist \
         -k $known_fusions \
         -o ${sample_id}.fusions.tsv \
         -O ${sample_id}.fusions.discarded.tsv
    """
}