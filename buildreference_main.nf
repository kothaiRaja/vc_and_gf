nextflow.enable.dsl = 2

process DOWNLOAD_REF_GENOME {
	tag "Download Reference genome"
    output:
    path "genome.fa"

    publishDir "/home/kothai/cq-git-sample/vc_and_gf/data/actual/reference", mode: 'copy'

    script:
    """
    if [[ -n "${params.genome_path}" && -f "${params.genome_path}" ]]; then
        echo "Using existing genome: ${params.genome_path}"
        ln -s "${params.genome_path}" genome.fa
    elif [[ -f "${params.ref_genome}" ]]; then
        echo "Previously downloaded genome found at ${params.ref_genome}"
        ln -s "${params.ref_genome}" genome.fa
    else
        echo "Downloading genome from ${params.genome_download_url}..."
        mkdir -p \$(dirname genome.fa)
        wget -O genome.fa.gz "${params.genome_download_url}"
        gunzip genome.fa.gz
    fi
    """
}

process DOWNLOAD_VARIANTS_SNP {
    tag "Download variants_snp VCF"

    output:
    path "variants_snp.vcf.gz"

    publishDir "/home/kothai/cq-git-sample/vc_and_gf/data/actual/reference", mode: 'copy'

    script:
    """
    if [[ -n "${params.variants_snp_path}" && -f "${params.variants_snp_path}" ]]; then
        echo " Using existing variants VCF: ${params.variants_snp_path}"
        ln -s "${params.variants_snp_path}" variants_snp.vcf.gz
    elif [[ -f "${params.ref_variants_snp}" ]]; then
        echo " Previously downloaded variants VCF found at ${params.ref_variants_snp}"
        ln -s "${params.ref_variants_snp}" variants_snp.vcf.gz
    else
        echo " Downloading variants VCF from ${params.variants_snp_download_url}..."
        mkdir -p \$(dirname variants_snp.vcf.gz)
        wget -O variants_snp.vcf.gz "${params.variants_snp_download_url}"
        echo " Download complete. Storing it permanently."
        mv variants_snp.vcf.gz "${params.ref_variants_snp}"
        ln -s "${params.ref_variants_snp}" variants_snp.vcf.gz
    fi
    """
}









workflow {

    // Step 1: Prepare the genome (either use existing or download)
    genome = DOWNLOAD_REF_GENOME()
	
	// Step 2: Download or use existing variants VCF
    variants_snp_vcf = DOWNLOAD_VARIANTS_SNP()

   // Print the output file paths
    println "✅ Reference genome is ready: ${genome}"
    println "✅ Variants VCF is ready: ${variants_snp_vcf}"
	
}