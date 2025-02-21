// ========================== Download Indels Variants VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_INDELS {
    tag "Check or Download Indels Variants"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz", emit: variants_indels

    when:
	!file("${params.test_data_dir}/reference/variants_indels.vcf.gz").exists()


    script:
    """
    wget -q -O variants_indels.vcf.gz ${params.variants_indels_download_url}
    """
}

process INDEX_INDEL_VCF {
    tag "Index INDEL VCF"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path vcf_file

    output:
    path "${vcf_file}.tbi", emit: indels_index

    script:
    """
    echo "Indexing INDEL VCF file: ${vcf_file}"
    tabix -p vcf ${vcf_file}
    """
}

// ========================== Download Indels Variants Index ========================== //
process DOWNLOAD_VARIANTS_INDELS_INDEX {
    tag "Download INDELS Index"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz.tbi", emit: indels_index

    script:
    """
    wget -q -O variants_indels.vcf.gz.tbi ${params.variants_indels_index_download_url}
    """
}