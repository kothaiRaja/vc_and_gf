nextflow.enable.dsl = 2


// ========================== Download Reference Genome ========================== //
process CHECK_OR_DOWNLOAD_REF_GENOME {
    tag "Check or Download Reference Genome"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa", emit: genome

    when:
	!file("${params.actual_data_dir}/reference/genome.fa").exists()
   
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

// ========================== Download SNP Variants VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_SNP {
    tag "Check or Download SNP Variants"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz", emit: variants_snp

    when:
    !file("${params.actual_data_dir}/reference/variants_snp.vcf.gz").exists()

    script:
    """
    wget -q -O variants_snp.vcf.gz ${params.variants_snp_download_url}
    """
}

// ========================== Download Indels Variants VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_INDELS {
    tag "Check or Download Indels Variants"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz", emit: variants_indels

    when:
	!file("${params.actual_data_dir}/reference/variants_indels.vcf.gz").exists()


    script:
    """
    wget -q -O variants_indels.vcf.gz ${params.variants_indels_download_url}
    """
}

// ========================== Download GTF Annotation File ========================== //
process CHECK_OR_DOWNLOAD_GTF {
    tag "Check or Download GTF File"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "annotations.gtf", emit: gtf

    script:
    """
    wget -q -O annotations.gtf ${params.gtf_download_url}
    """
}

// ========================== Download Denylist File ========================== //
process CHECK_OR_DOWNLOAD_DENYLIST {
    tag "Check or Download Denylist BED File"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed", emit: denylist

    script:
    """
    wget -q -O denylist.bed ${params.denylist_download_url}
    """
}

// ========================== Create FASTA Index ========================== //
process CREATE_FASTA_INDEX {
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input:
    path genome_fasta

    output:
    path "${genome_fasta}.fai", emit: fasta_index

    script:
    """
    samtools faidx ${genome_fasta}
    """
}

// ========================== Create Genome Dictionary ========================== //
process CREATE_GENOME_DICT {
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input:
    path genome_fasta

    output:
    path "${genome_fasta.baseName}.dict", emit: genome_dict

    script:
    """
    gatk CreateSequenceDictionary -R $genome_fasta -O ${genome_fasta.baseName}.dict
    """
}

// ========================== Create STAR Genome Index ========================== //
process CREATE_STAR_INDEX {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

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

process PREPARE_VCF_FILE {
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input: 
    path snpsFile
    path indelsFile
    path denylisted

    output:
    tuple path("merged.filtered.recode.vcf.gz"), 
          path("merged.filtered.recode.vcf.gz.tbi"), emit: filtered_vcf


		  
	
    script:  
    """
    # Filter SNPs file
    bcftools view -T ^${denylisted} ${snpsFile} -Oz -o ${snpsFile.baseName}.filtered.recode.vcf.gz
    tabix -p vcf ${snpsFile.baseName}.filtered.recode.vcf.gz

    # Filter INDELs file
    bcftools view -T ^${denylisted} ${indelsFile} -Oz -o ${indelsFile.baseName}.filtered.recode.vcf.gz
    tabix -p vcf ${indelsFile.baseName}.filtered.recode.vcf.gz

    # Merge the filtered SNPs and INDELs into a single VCF file
    bcftools merge ${snpsFile.baseName}.filtered.recode.vcf.gz ${indelsFile.baseName}.filtered.recode.vcf.gz \
        -Oz -o merged.filtered.recode.vcf.gz

    # Create a tabix index for the merged VCF
    tabix -p vcf merged.filtered.recode.vcf.gz
	
	# Check and fix contig prefixes in the merged VCF file
    if zcat merged.filtered.recode.vcf.gz | grep -q "^##contig=<ID=chr"; then
        echo "Renaming contig prefixes..."
        zcat merged.filtered.recode.vcf.gz | sed 's/^##contig=<ID=chr/##contig=<ID=/' | bgzip > fixed_merged.filtered.recode.vcf.gz
        tabix -p vcf fixed_merged.filtered.recode.vcf.gz
        mv fixed_merged.filtered.recode.vcf.gz merged.filtered.recode.vcf.gz
        mv fixed_merged.filtered.recode.vcf.gz.tbi merged.filtered.recode.vcf.gz.tbi
    else
        echo "Contig prefixes are already correct."
    fi
    """
}

process CHECK_JAVA {
    tag "Check Java"
    container null
	
    output:
    path "java_check.log", emit: java_output
	
    when:
    true

    script:
    """
    # Get the Java version
    java_version=\$(java -version 2>&1 | grep -oP '(?<=version ")([0-9]+)' | head -1)
    
    if [[ -z "\$java_version" ]]; then
        echo "ERROR: Java is not installed or not in PATH." > java_check.log
    elif [[ \$java_version -lt 21 ]]; then
        echo "WARNING: Java version \$java_version detected. Please update to version 21 or higher." > java_check.log
    else
        echo "Java version \$java_version detected. Java is up-to-date." > java_check.log
    fi
    """
}




process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
	publishDir "${params.actual_data_dir}/Tools", mode: 'copy'
	container null
	

	  
    output:
    path "${params.snpeff_jar_dir}/snpEff.jar"
	path "${params.snpeff_jar_dir}/snpEff.config"
	
	script:
    """
    mkdir -p ${params.snpeff_jar_dir}
    wget -q -O snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -j snpEff_latest_core.zip -d ${params.snpeff_jar_dir}
    rm snpEff_latest_core.zip
	
    """
}

process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
	publishDir "${params.actual_data_dir}/Tools/snpEff", mode: 'copy'
	container null
	
	
    input:
    val genome
    path snpeff_jar_path

    output:
    path "${params.snpeff_db_dir}/${genome}"
	
	script:
    """
    # Ensure the output directory exists first
    mkdir -p ${params.snpeff_db_dir}

    # Use an absolute path for the data directory
    data_dir=\$(realpath ${params.snpeff_db_dir})

    # Download the database
    # Check if the database is already downloaded
    if [ ! -d "\$data_dir/${genome}" ]; then
        echo "Downloading SnpEff database for ${genome}..."
        java -Xmx4g -Xms2g -jar ${snpeff_jar_path} download ${genome} -dataDir \$data_dir -v
    else
        echo "SnpEff database for ${genome} already exists. Skipping download."
    fi
    """
}

process DOWNLOAD_ARRIBA {
    container null
	publishDir "${params.actual_data_dir}/Tools/ARRIBA", mode: 'copy' 
    
	output:
    path 'arriba_v2.4.0', emit: 'arriba_dir'
	
	script:
    """
    # Define the Arriba download URL and target directory
    URL="https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz"
    TARGET_DIR="arriba_v2.4.0"

    # Download Arriba tarball
    wget -O arriba.tar.gz \$URL

    # Create the target directory and extract the tarball
    mkdir -p \$TARGET_DIR
    tar -xzvf arriba.tar.gz -C \$TARGET_DIR --strip-components=1

    # Clean up the downloaded tarball
    rm arriba.tar.gz

    # Output the extracted directory
    echo "Arriba tool extracted to: \$TARGET_DIR"
    """
}

process DOWNLOAD_VEP_CACHE {
    tag "Downloading VEP Cache"
	container null
	publishDir "${params.actual_data_dir}/Tools/VEP", mode: 'copy'
	
	output:
    path "${params.vep_cache}/homo_sapiens/110_GRCh38", type: 'dir'


    script:
    """
     mkdir -p vep_cache
    wget -O homo_sapiens_vep_110_GRCh38.tar.gz https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz
    
    # Check if the download was successful
    if [ ! -s homo_sapiens_vep_110_GRCh38.tar.gz ]; then
        echo "Error: Failed to download VEP cache file" >&2
        exit 1
    fi

    tar -xvzf homo_sapiens_vep_110_GRCh38.tar.gz -C vep_cache
    rm homo_sapiens_vep_110_GRCh38.tar.gz
    """
}

process DOWNLOAD_CLINVAR {
    tag "Downloading ClinVar VCF"
	container null
	publishDir "${params.actual_data_dir}/Tools/VEP", mode: 'copy'
	
	 when:
    !file("${params.clinvar}").exists() || !file("${params.clinvartbi}").exists()
	
    output:
    path "clinvar.vcf.gz"
    path "clinvar.vcf.gz.tbi"

    script:
    """
    wget -O clinvar.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
    
    # Check if the VCF download was successful
    if [ ! -s clinvar.vcf.gz ]; then
        echo "Error: Failed to download ClinVar VCF file" >&2
        exit 1
    fi

    wget -O clinvar.vcf.gz.tbi ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
    
    # Check if the TBI index download was successful
    if [ ! -s clinvar.vcf.gz.tbi ]; then
        echo "Error: Failed to download ClinVar VCF index file" >&2
        exit 1
    fi
    """
}

process WRITE_REFERENCE_CONFIG {
    tag "Generate reference_paths.config"
    publishDir "${params.base_dir}", mode: 'copy'

    input:
        path genome
        path variants_snp 
        path variants_indels 
        path gtf 
        path denylist
        path fasta_index 
        path genome_dict 
        path star_index 
        tuple path(filtered_vcf), path(filtered_vcf_tbi)

    output:
        path "reference_paths.config"

    script:
    """
    echo "🔹 Debugging Input Paths:"
    echo "Genome: ${genome}"
    echo "SNP Variants: ${variants_snp}"
    echo "INDEL Variants: ${variants_indels}"
    echo "GTF: ${gtf}"
    echo "Denylist: ${denylist}"
    echo "FASTA Index: ${fasta_index}"
    echo "Genome Dict: ${genome_dict}"
    echo "STAR Index: ${star_index}"
    echo "Filtered VCF: ${filtered_vcf}"
    echo "Filtered VCF TBI: ${filtered_vcf_tbi}"

    # Ensure the reference directory exists
    mkdir -p ${params.base_dir}

    # Write reference paths config
    cat <<EOL > reference_paths.config
	params {
    ref_genome = "\$(realpath ${params.actual_data_dir}/reference/genome.fa)"
    ref_variants_snp = "\$(realpath ${params.actual_data_dir}/reference/variants_snp.vcf.gz)"
    ref_variants_indels = "\$(realpath ${params.actual_data_dir}/reference/variants_indels.vcf.gz)"
    ref_gtf = "\$(realpath ${params.actual_data_dir}/reference/annotations.gtf)"
    ref_denylist = "\$(realpath ${params.actual_data_dir}/reference/denylist.bed)"
    ref_fasta_index = "\$(realpath ${params.actual_data_dir}/reference/genome.fa.fai)"
    ref_genome_dict = "\$(realpath ${params.actual_data_dir}/reference/genome.dict)"
    ref_star_index = "\$(realpath ${params.actual_data_dir}/reference/STAR_index)"
    ref_filtered_vcf = "\$(realpath ${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz)"
    ref_filtered_vcf_tbi = "\$(realpath ${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi)"
}
EOL
	
	# Verify that the file was created
    if [ ! -f ${params.base_dir}/reference_paths.config ]; then
        echo "❌ ERROR: reference_paths.config was not created!"
        exit 1
    else
        echo "✅ reference_paths.config successfully created!"
    fi

    # List directory to check if file was created
    ls -lh ${params.actual_data_dir}/reference
    """
}








// ========================== Helper Function ========================== //
def safeFileChannel(localPath, centralPath, processFunc) {
    if (file(localPath).exists()) {
        log.info "✅ Using existing file: ${localPath}"
        return Channel.value(file(localPath))  
    } 
    else if (centralPath && file(centralPath).exists()) {
        log.info "✅ Using centralized file: ${centralPath}"
        return Channel.value(file(centralPath))  
    } 
    else {
        log.info "📥 File not found, running process to create it..."
        return processFunc().map { it -> 
            def publishedPath = file("${params.actual_data_dir}/reference/${it.getName()}") 
            log.info "✅ Published file path: ${publishedPath}"
            return publishedPath
        }
    }
}


// ========================== Workflow Definition ========================== //
workflow {
    log.info "🚀 Starting the pipeline with reference files processing workflow..."

    def genome_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/genome.fa",
        params.central_genome_path,
        CHECK_OR_DOWNLOAD_REF_GENOME
    )

    def variants_snp_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/variants_snp.vcf.gz",
        params.central_variants_snp_path,
        CHECK_OR_DOWNLOAD_VARIANTS_SNP
    )

    def variants_indels_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/variants_indels.vcf.gz",
        params.central_variants_indels_path,
        CHECK_OR_DOWNLOAD_VARIANTS_INDELS
    )

    def gtf_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/annotations.gtf",
        params.central_gtf_path,
        CHECK_OR_DOWNLOAD_GTF
    )

    def denylist_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/denylist.bed",
        params.central_denylist_path,
        CHECK_OR_DOWNLOAD_DENYLIST
    )

    def fasta_index_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/genome.fa.fai",
        params.central_fasta_index_path,
        { CREATE_FASTA_INDEX(genome_channel) }
    )

    def genome_dict_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/genome.dict",
        params.central_genome_dict_path,
        { CREATE_GENOME_DICT(genome_channel) }
    )

    def star_index_channel = safeFileChannel(
        "${params.actual_data_dir}/reference/STAR_index",
        params.central_star_index_path,
        { CREATE_STAR_INDEX(genome_channel, gtf_channel) }
    )

	def filtered_vcf_channel = PREPARE_VCF_FILE(
        variants_snp_channel, 
        variants_indels_channel, 
        denylist_channel
    )

    // Extract tuple values correctly
    def filtered_vcf_path = filtered_vcf_channel.map { file("${params.actual_data_dir}/reference/" + it[0].getName()) }
    def filtered_vcf_tbi_path = filtered_vcf_channel.map { file("${params.actual_data_dir}/reference/" + it[1].getName()) }




    // ========================== ✅ WRITE_REFERENCE_CONFIG Call ========================== //
    WRITE_REFERENCE_CONFIG(
			genome_channel,
            variants_snp_channel,
            variants_indels_channel,
            gtf_channel,
            denylist_channel,
            fasta_index_channel,
            genome_dict_channel,
            star_index_channel,
			filtered_vcf_channel
            
        
    )

    log.info "🎉 Reference files processing workflow completed!"
}




 
    
	

   


