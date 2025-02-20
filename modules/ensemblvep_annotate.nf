process ANNOTATE_INDIVIDUAL_VARIANTS_VEP {
    tag "${sample_id}_vep_annotate"

    container "https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index), val(strandedness) // Input VCF and its index
    path vep_cache                                                  // VEP cache directory
    path clinvar_vcf                                                // ClinVar VCF file
    path clinvar_index                                              // ClinVar Tabix index file

    output:
    tuple val(sample_id), path("${sample_id}.vep.annotated.vcf"),val(strandedness), path("${sample_id}.vep.summary.html")

    script:
    """
    # Annotate using Ensembl VEP with ClinVar
    vep --input_file ${filtered_vcf} \
        --output_file ${sample_id}.vep.annotated.vcf \
        --stats_file ${sample_id}.vep.summary.html \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN

    # Validate that the annotated VCF is not empty
    if [ ! -s ${sample_id}.vep.annotated.vcf ]; then
        echo "Error: VEP annotation output is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}

process ANNOTATEVARIANTS_VEP {
    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path input_vcf          // Input VCF file
    path input_vcf_tbi      // Tabix index file for the VCF
	path tsv
    path vep_cache          // VEP cache directory
    path clinvar_vcf        // ClinVar VCF file
    path clinvar_vcf_tbi    // Tabix index file for ClinVar
	

    output:
    path "annotated_variants.vcf", emit: annotated_vcf  // ✅ Ensures correct emission!
    path "annotated_variants.html", emit: summary_html

    script:
    """
    # Debugging: Log resolved paths and existence
    echo "Input VCF: \$(realpath ${input_vcf})"
    echo "Input VCF Index: \$(realpath ${input_vcf_tbi})"
    echo "VEP Cache Directory: \$(realpath ${vep_cache})"
    echo "ClinVar VCF: \$(realpath ${clinvar_vcf})"
    echo "ClinVar Index: \$(realpath ${clinvar_vcf_tbi})"

    if [ ! -d "${vep_cache}" ]; then
        echo "Error: VEP cache directory does not exist at: \$(realpath ${vep_cache})" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf}" ]; then
        echo "Error: ClinVar VCF file not found at: \$(realpath ${clinvar_vcf})" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf_tbi}" ]; then
        echo "Error: ClinVar index file not found at: \$(realpath ${clinvar_vcf_tbi})" >&2
        exit 1
    fi

    # Run VEP annotation
    vep \
        --input_file ${input_vcf} \
        --output_file annotated_variants.vcf \
        --stats_file annotated_variants.html \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN
    """
}