nextflow.enable.dsl = 2


process DOWNLOAD_REF_GENOME {
    tag "Download test genome"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa"

    when:
    !file("${params.test_data_dir}/reference/genome.fa").exists()

    script:
    """
    wget -q -O genome.fa ${params.test_data_genome}
    """
}

process DOWNLOAD_VARIANTS_SNP {
    tag "Download test variants VCF"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz"

    when:
    !file("${params.test_data_dir}/reference/variants_snp.vcf.gz").exists()

    script:
    """
    wget -q -O variants_snp.vcf.gz ${params.test_data_dbsnp}
    """
}

process DOWNLOAD_VARIANTS_INDELS {
    tag "Download test variants VCF"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz"

    when:
    !file("${params.test_data_dir}/reference/variants_indels.vcf.gz").exists()

    script:
    """
    wget -q -O variants_indels.vcf.gz ${params.test_data_known_indels}
    """
}

process DOWNLOAD_DENYLIST {
    tag "Download test denylist BED"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed"

    when:
    !file("${params.test_data_dir}/reference/denylist.bed").exists()

    script:
    """
    wget -q -O denylist.bed ${params.test_data_denylist}
    """
}

process DOWNLOAD_GTF {
    tag "Download test GTF"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "annotations.gtf"

    when:
    !file("${params.test_data_dir}/reference/annotations.gtf").exists()

    script:
    """
    wget -q -O annotations.gtf ${params.test_data_gtf}
    """
}

//Create a fasta Genome Index

process CREATE_FASTA_INDEX {
	container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'
	
    input:
    path genome_fasta

    output:
    path("${genome_fasta}.fai")

    script:
    """
    samtools faidx ${genome_fasta}
    """
}

// Process: Create Genome Dictionary

process CREATE_GENOME_DICT {
	container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'
    input:
    path genome_fasta

    output:
    path("${genome_fasta.baseName}.dict")

    script:
    """
    gatk CreateSequenceDictionary -R $genome_fasta -O ${genome_fasta.baseName}.dict
    """
}

//Process: STAR Genome Index

process CREATE_STAR_INDEX {
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path genome_fasta
    path genome_gtf

    output:
    path "STAR_index", type: 'dir'

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
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input: 
    path snpsFile
    path indelsFile
    path denylisted

    output:
    tuple path("merged.filtered.recode.vcf.gz"),
          path("merged.filtered.recode.vcf.gz.tbi")

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
    stdout emit: java_output

    script:
    """
    # Get the Java version
    java_version=\$(java -version 2>&1 | grep -oP '(?<=version ")([0-9]+)' | head -1)
    
    if [[ -z "\$java_version" ]]; then
        echo "ERROR: Java is not installed or not in PATH."
        exit 1
    elif [[ \$java_version -lt 21 ]]; then
        echo "WARNING: Java version \$java_version detected. Please update to version 21 or higher."
        exit 1
    else
        echo "Java version \$java_version detected. Java is up-to-date."
    fi
    """
}



process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
	publishDir "${params.test_data_dir}/Tools", mode: 'copy'
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
	publishDir "${params.test_data_dir}/Tools/snpEff", mode: 'copy'
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
	publishDir "${params.test_data_dir}/Tools/ARRIBA", mode: 'copy' 
    
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
	publishDir "${params.test_data_dir}/Tools/VEP", mode: 'copy'
	
	when:
    !file("${params.vep_cache}/homo_sapiens/110_GRCh38").exists()

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
	publishDir "${params.test_data_dir}/Tools/VEP", mode: 'copy'
	
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




workflow {
    // Step 1: Check and download reference genome
    def genome_path = "${params.test_data_dir}/reference/genome.fa"
    if (file(genome_path).exists()) {
        log.info "Reference genome already exists at: ${genome_path}"
    } else {
        log.info "Downloading reference genome..."
        DOWNLOAD_REF_GENOME()
    }

    // Step 2: Check and download SNPs VCF file
    def snps_path = "${params.test_data_dir}/reference/variants_snp.vcf.gz"
    if (file(snps_path).exists()) {
        log.info "SNPs VCF file already exists at: ${snps_path}"
    } else {
        log.info "Downloading SNPs VCF file..."
        DOWNLOAD_VARIANTS_SNP()
    }

    // Step 3: Check and download INDELs VCF file
    def indels_path = "${params.test_data_dir}/reference/variants_indels.vcf.gz"
    if (file(indels_path).exists()) {
        log.info "INDELs VCF file already exists at: ${indels_path}"
    } else {
        log.info "Downloading INDELs VCF file..."
        DOWNLOAD_VARIANTS_INDELS()
    }

    // Step 4: Check and download denylist BED file
    def denylist_path = "${params.test_data_dir}/reference/denylist.bed"
    if (file(denylist_path).exists()) {
        log.info "Denylist BED file already exists at: ${denylist_path}"
    } else {
        log.info "Downloading denylist BED file..."
        DOWNLOAD_DENYLIST()
    }

    // Step 5: Check and download GTF file
    def gtf_path = "${params.test_data_dir}/reference/annotations.gtf"
    if (file(gtf_path).exists()) {
        log.info "GTF file already exists at: ${gtf_path}"
    } else {
        log.info "Downloading GTF file..."
        DOWNLOAD_GTF()
    }

    // Step 6: Check and create genome fasta index
    def fasta_index_path = "${params.test_data_dir}/reference/genome.fa.fai"
    if (file(fasta_index_path).exists()) {
        log.info "Fasta index already exists at: ${fasta_index_path}"
    } else {
        log.info "Creating fasta index..."
        CREATE_FASTA_INDEX(file(genome_path))
    }

    // Step 7: Check and create genome dictionary
    def genome_dict_path = "${params.test_data_dir}/reference/genome.dict"
    if (file(genome_dict_path).exists()) {
        log.info "Genome dictionary already exists at: ${genome_dict_path}"
    } else {
        log.info "Creating genome dictionary..."
        CREATE_GENOME_DICT(file(genome_path))
    }

    // Step 8: Check and create STAR index
    def star_index_path = "${params.test_data_dir}/reference/STAR_index"
    if (file(star_index_path).exists()) {
        log.info "STAR index already exists at: ${star_index_path}"
    } else {
        log.info "Creating STAR index..."
        CREATE_STAR_INDEX(file(genome_path), file(gtf_path))
    }

    // Step 9: Check and prepare VCF file (merge SNPs and INDELs)
    def filtered_vcf_path = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz"
    def filtered_vcf_index_path = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"
    if (file(filtered_vcf_path).exists() && file(filtered_vcf_index_path).exists()) {
        log.info "Filtered VCF file already exists at: ${filtered_vcf_path}"
    } else {
        log.info "Preparing filtered VCF file (merging SNPs and INDELs)..."
        PREPARE_VCF_FILE(file(snps_path), file(indels_path), file(denylist_path))
    }

    // Step 10: Check Java installation (always executed)
    log.info "Checking Java installation..."
    // Call the CHECK_JAVA process
    def java_check_result = CHECK_JAVA()

    // Print the Java version check result to the console
    java_check_result.java_output.view()

    //log.info "Checking SnpEff tool existence..."
	def snpeff_jar_path = "${params.test_data_dir}/Tools/snpEff/snpEff.jar"
	def snpeff_config_path = "${params.test_data_dir}/Tools/snpEff/snpEff.config"

	// Check if SnpEff tool files exist
	if (file(snpeff_jar_path).exists() && file(snpeff_config_path).exists()) {
    log.info "SnpEff tool already exists:"
    log.info "  - JAR file: ${snpeff_jar_path}"
    log.info "  - Config file: ${snpeff_config_path}"
	} else {
    log.info "SnpEff tool not found. Downloading..."
    DOWNLOAD_SNPEFF_TOOL()
	}


	log.info "Checking SnpEff database existence..."
	def snpeff_db_path = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"

	// Check if the SnpEff database exists
	if (file(snpeff_db_path).exists()) {
    log.info "SnpEff database already exists at: ${snpeff_db_path}"
	} else {
    log.info "SnpEff database not found for genome version '${params.genomedb}'. Downloading..."
    DOWNLOAD_SNPEFF_DB(params.genomedb, file(snpeff_jar_path))
	}


    // Step 13: Check and download Arriba tool
    def arriba_dir_path = "${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0"
    if (file(arriba_dir_path).exists()) {
        log.info "Arriba tool already exists at: ${arriba_dir_path}"
    } else {
        log.info "Downloading Arriba tool..."
        DOWNLOAD_ARRIBA()
    }

    // Step 14: Check and download VEP cache
    def vep_cache_path = "${params.test_data_dir}/Tools/VEP/vep_cache"
    if (file(vep_cache_path).exists()) {
        log.info "VEP cache already exists at: ${vep_cache_path}"
    } else {
        log.info "Downloading VEP cache..."
        DOWNLOAD_VEP_CACHE()
    }

    // Step 15: Check and download ClinVar files
    def clinvar_vcf_path = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz"
    def clinvar_tbi_path = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi"
    if (file(clinvar_vcf_path).exists() && file(clinvar_tbi_path).exists()) {
        log.info "ClinVar VCF files already exist at: ${clinvar_vcf_path} and ${clinvar_tbi_path}"
    } else {
        log.info "Downloading ClinVar VCF files..."
        DOWNLOAD_CLINVAR()
    }

    // Log workflow completion
    log.info "Reference files processing workflow completed successfully!"
}


