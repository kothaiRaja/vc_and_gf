// ========================== Create STAR Genome Index ========================== //
process CREATE_STAR_INDEX {
    tag "Create STAR Genome Index"
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path genome_fasta
    path genome_gtf

    output:
    path "STAR_index", type: 'dir', emit: STAR_index

    script:
    """
    mkdir -p STAR_index
    STAR --runMode genomeGenerate \
         --genomeDir STAR_index \
         --genomeFastaFiles ${genome_fasta} \
         --sjdbGTFfile ${genome_gtf} \
         --runThreadN 8
    """
}