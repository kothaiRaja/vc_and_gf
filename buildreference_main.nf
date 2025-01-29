nextflow.enable.dsl = 2

process DOWNLOAD_REF_GENOME {
    tag "Download Reference Genome"

    output:
    path "genome.fa"

    publishDir "/home/kothai/cq-git-sample/vc_and_gf/data/actual/reference", mode: 'copy'

    script:
    """
    if [[ -n "${params.genome_path}" && -f "${params.genome_path}" ]]; then
        echo " Using existing reference genome: ${params.genome_path}"
        ln -s "${params.genome_path}" genome.fa
    elif [[ -f "${params.ref_genome}" ]]; then
        echo " Previously downloaded reference genome found at ${params.ref_genome}"
        ln -s "${params.ref_genome}" genome.fa
    else
        echo " Downloading reference genome from ${params.genome_download_url}..."
        mkdir -p "\$(dirname ${params.ref_genome})"
        wget -O genome.fa.gz "${params.genome_download_url}"
        gunzip genome.fa.gz
        echo " Download complete. Storing it permanently at ${params.ref_genome}"
        mv genome.fa "${params.ref_genome}"
        ln -s "${params.ref_genome}" genome.fa
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

process DOWNLOAD_VARIANTS_INDELS {
    tag "Download variants Indels VCF"

    output:
    path "variants_indels.vcf.gz"

    publishDir "/home/kothai/cq-git-sample/vc_and_gf/data/actual/reference", mode: 'copy'

    script:
    """
    if [[ -n "${params.variants_indels_path}" && -f "${params.variants_indels_path}" ]]; then
        echo " Using existing variants Indels VCF: ${params.variants_indels_path}"
        ln -s "${params.variants_indels_path}" variants_indels.vcf.gz
    elif [[ -f "${params.ref_variants_indels}" ]]; then
        echo " Previously downloaded variants Indels VCF found at ${params.ref_variants_indels}"
        ln -s "${params.ref_variants_indels}" variants_indels.vcf.gz
    else
        echo " Downloading variants Indels VCF from ${params.variants_indels_download_url}..."
        mkdir -p "\$(dirname ${params.ref_variants_indels})"
        wget -O variants_indels.vcf.gz "${params.variants_indels_download_url}"
        echo " Download complete. Storing it permanently at ${params.ref_variants_indels}"
        mv variants_indels.vcf.gz "${params.ref_variants_indels}"
        ln -s "${params.ref_variants_indels}" variants_indels.vcf.gz
    fi
    """
}

process DOWNLOAD_GTF {
    tag "Download GTF annotation file"

    output:
    path "annotations.gtf"

    publishDir "/home/kothai/cq-git-sample/vc_and_gf/data/actual/reference", mode: 'copy'

    script:
    """
    if [[ -n "${params.gtf_path}" && -f "${params.gtf_path}" ]]; then
        echo " Using existing GTF file: ${params.gtf_path}"
        ln -s "${params.gtf_path}" annotations.gtf
    elif [[ -f "${params.ref_gtf}" ]]; then
        echo " Previously downloaded GTF file found at ${params.ref_gtf}"
        ln -s "${params.ref_gtf}" annotation.gtf
    else
        echo " Downloading GTF file from ${params.gtf_download_url}..."
        mkdir -p "\$(dirname ${params.ref_gtf})"
        wget -O annotations.gtf.gz "${params.gtf_download_url}"
        echo " Download complete. Unzipping..."
        gunzip -c annotations.gtf.gz > annotations.gtf
        echo " Unzipped and storing permanently at ${params.ref_gtf}"
        mv annotations.gtf "${params.ref_gtf}"
        ln -s "${params.ref_gtf}" annotations.gtf
    fi
    """
}


process DOWNLOAD_DENY_LIST {
    tag "Download Deny List BED file"

    output:
    path "denylist.bed"

    publishDir "/home/kothai/cq-git-sample/vc_and_gf/data/actual/reference", mode: 'copy'

    script:
    """
    if [[ -n "${params.denylist_path}" && -f "${params.denylist_path}" ]]; then
        echo "✅ Using existing Deny List BED file: ${params.denylist_path}"
        ln -s "${params.denylist_path}" denylist.bed
    elif [[ -f "${params.ref_denylist}" ]]; then
        echo "🔄 Previously downloaded Deny List file found at ${params.ref_denylist}"
        ln -s "${params.ref_denylist}" denylist.bed
    else
        echo "🚀 Downloading Deny List BED file from ${params.denylist_download_url}..."
        mkdir -p "\$(dirname ${params.ref_denylist})"
        wget -O denylist.bed.gz "${params.denylist_download_url}"
        echo "✅ Download complete. Unzipping..."
        gunzip -c denylist.bed.gz > denylist.bed
        echo "✅ Unzipped and storing permanently at ${params.ref_denylist}"
        mv denylist.bed "${params.ref_denylist}"
        ln -s "${params.ref_denylist}" denylist.bed
    fi
    """
}









workflow {

    // Step 1: Prepare the genome (either use existing or download)
    genome = DOWNLOAD_REF_GENOME()
	
	// Step 2: Download or use existing variants VCF
    variants_snp = DOWNLOAD_VARIANTS_SNP()


    // Step 3: Download variants Indels VCF
    variants_indels = DOWNLOAD_VARIANTS_INDELS()

    // Step 4: Download GTF annotation file
    gtf_file = DOWNLOAD_GTF()

    // Step 5: Download Deny List BED file
    deny_list = DOWNLOAD_DENY_LIST()

    // Print outputs
    println "✅ Reference genome ready: ${genome}"
    println "✅ Variants SNP VCF ready: ${variants_snp}"
    println "✅ Variants Indels VCF ready: ${variants_indels}"
    println "✅ GTF annotation file ready: ${gtf_file}"
    println "✅ Deny List BED file ready: ${deny_list}"
	
}