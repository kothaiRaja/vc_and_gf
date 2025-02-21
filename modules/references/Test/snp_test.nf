// ========================== Download SNP Variants VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_SNP {
    tag "Check or Download SNP Variants"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz", emit: variants_snp

    when:
    !file("${params.test_data_dir}/reference/variants_snp.vcf.gz").exists()

    script:
    """
    wget -q -O variants_snp.vcf.gz ${params.variants_snp_download_url}
    """
}

// ========================== Download SNP Variants Index ========================== //
process DOWNLOAD_VARIANTS_SNP_INDEX {
    tag "Download SNP Index"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz.tbi", emit: snp_index

    script:
    """
    wget -q -O variants_snp.vcf.gz.tbi ${params.variants_snp_index_download_url}
    """
}

//==========================Create Indexing=============================//
process INDEX_SNP_VCF {
    tag "Index SNP VCF"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path vcf_file

    output:
    path "${vcf_file}.tbi", emit: snp_index

    script:
    """
    echo "Indexing SNP VCF file: ${vcf_file}"
    tabix -p vcf ${vcf_file}
    """
}