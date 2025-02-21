// ========================== Download Reference Genome and Index ========================== //
process CHECK_OR_DOWNLOAD_REF_GENOME {
    tag "Check or Download Reference Genome"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa", emit: genome

    when:
	!file("${params.test_data_dir}/reference/genome.fa").exists()
   
	script:
    """
    wget -q -O genome.fa.gz ${params.genome_download_url}

    if file genome.fa.gz | grep -q 'gzip'; then
        gunzip genome.fa.gz
    else
        mv genome.fa.gz genome.fa
    fi
    """
}