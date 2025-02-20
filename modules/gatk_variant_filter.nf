process GATK_VARIANT_FILTER {
    tag "${sample_id}_variant_filter"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variant_filter/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(strandedness)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), val(strandedness)


    script:
    """
    # Run GATK VariantFiltration with configurable thresholds
    gatk VariantFiltration \
        -R ${genome} \
        -V ${vcf_file} \
        --cluster-window-size ${params.gatk_vf_window_size} \
        --cluster-size ${params.gatk_vf_cluster_size} \
        --filter-name "LowQual" --filter-expression "QUAL < ${params.gatk_vf_qual_filter}" \
        --filter-name "LowQD" --filter-expression "QD < ${params.gatk_vf_qd_filter}" \
        --filter-name "HighFS" --filter-expression "FS > ${params.gatk_vf_fs_filter}" \
        --filter-name "LowMQ" --filter-expression "MQ < ${params.gatk_vf_mq_filter}" \
        --filter-name "HighSOR" --filter-expression "SOR > ${params.gatk_vf_sor_filter}" \
        --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < ${params.gatk_vf_read_pos_filter}" \
        --filter-name "LowBaseQRankSum" --filter-expression "BaseQRankSum < ${params.gatk_vf_baseq_filter}" \
        -O ${sample_id}_filtered.vcf.gz

    # Index the filtered VCF file for downstream compatibility
    gatk IndexFeatureFile -I ${sample_id}_filtered.vcf.gz

    # Redirect logs
    if [ ! -s ${sample_id}_filtered.vcf.gz ]; then
        echo "Error: Filtered VCF is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}