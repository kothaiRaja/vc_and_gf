//==================================Download Genome Index======================================//


process DOWNLOAD_GENOME_INDEX {
    tag "Download Genome Index"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa.fai", emit: genome_fai

    script:
    """
    echo "⚠️ Downloading genome index from provided URL..."
    wget -q -O genome.fa.fai ${params.genome_index_download_url}
    """
}

//====================================Create Genome Index===================================// 
process CREATE_GENOME_INDEX {
    tag "Create Genome Index"
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.14--hb421002_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa

    output:
    path "genome.fa.fai", emit: genome_fai

    script:
    """
    echo "⚠️ Creating genome index using samtools..."
    samtools faidx $genome_fa
    """
}