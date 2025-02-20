process BCFTOOLS_QUERY {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
    publishDir "${params.outdir}/variant_summary", mode: "copy"

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_index), val(strandedness)

    output:
    tuple val(sample_id), val(strandedness), path("filtered_variants_summary_${sample_id}.txt")

    script:
    """
    # Generate a summary of filtered variants
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${filtered_vcf} > filtered_variants_summary_${sample_id}.txt

   # Handle cases where the output file is empty
    if [ ! -s filtered_variants_summary_${sample_id}.txt ]; then
        echo "No variants detected for ${sample_id}" > filtered_variants_summary_${sample_id}.txt
        exit 0
    fi
    """
}