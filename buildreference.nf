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
   # Download the genome file
	wget -q -O genome.fa.gz ${params.test_data_genome}

	# Check if the file is gzip-compressed
	if file genome.fa.gz | grep -q 'gzip'; then
		gunzip genome.fa.gz
	else
		mv genome.fa.gz genome.fa
	fi
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
	
	when:
    !file("${genome_fasta}.fai").exists()

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
	
	when:
    !file("${genome_fasta.baseName}.dict").exists()

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
	
	 when:
    !file("${params.test_data_dir}/reference/STAR_index").exists()

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
		  
	when:
    !(file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz").exists() &&
      file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi").exists())

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
	
	when:
    true

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
	
	when:
    !(file("${params.snpeff_jar_dir}/snpEff.jar").exists() && file("${params.snpeff_jar_dir}/snpEff.config").exists())

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
	
	when:
    !file("${params.snpeff_db_dir}/${genome}").exists()

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
	
	when:
    !file("${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0").exists()

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
    log.info "Starting reference files processing workflow..."

    // Step 1: Reference genome
    def genome_path = "${params.test_data_dir}/reference/genome.fa"
    def genome_channel = file(genome_path).exists() 
        ? Channel.value(file(genome_path)) 
        : DOWNLOAD_REF_GENOME()

    // Step 2: SNPs VCF file
    def snps_path = "${params.test_data_dir}/reference/variants_snp.vcf.gz"
    def snps_channel = file(snps_path).exists() 
        ? Channel.value(file(snps_path)) 
        : DOWNLOAD_VARIANTS_SNP()

    // Step 3: INDELs VCF file
    def indels_path = "${params.test_data_dir}/reference/variants_indels.vcf.gz"
    def indels_channel = file(indels_path).exists() 
        ? Channel.value(file(indels_path)) 
        : DOWNLOAD_VARIANTS_INDELS()

    // Step 4: Denylist BED file
    def denylist_path = "${params.test_data_dir}/reference/denylist.bed"
    def denylist_channel = file(denylist_path).exists() 
        ? Channel.value(file(denylist_path)) 
        : DOWNLOAD_DENYLIST()

    // Step 5: GTF file
    def gtf_path = "${params.test_data_dir}/reference/annotations.gtf"
    def gtf_channel = file(gtf_path).exists() 
        ? Channel.value(file(gtf_path)) 
        : DOWNLOAD_GTF()

    // Step 6: Fasta index
    def fasta_index_path = "${params.test_data_dir}/reference/genome.fa.fai"
    def fasta_index_channel = file(fasta_index_path).exists() 
        ? Channel.value(file(fasta_index_path)) 
        : CREATE_FASTA_INDEX(genome_channel)

    // Step 7: Genome dictionary
    def genome_dict_path = "${params.test_data_dir}/reference/genome.dict"
    def genome_dict_channel = file(genome_dict_path).exists() 
        ? Channel.value(file(genome_dict_path)) 
        : CREATE_GENOME_DICT(genome_channel)

    // Step 8: STAR index
    def star_index_path = "${params.test_data_dir}/reference/STAR_index"
    def star_index_channel = file(star_index_path).exists() 
        ? Channel.value(file(star_index_path)) 
        : CREATE_STAR_INDEX(genome_channel, gtf_channel)

    // Step 9: Filtered VCF file
    def filtered_vcf_path = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz"
    def filtered_vcf_index_path = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"
    def filtered_vcf_channel = (file(filtered_vcf_path).exists() && file(filtered_vcf_index_path).exists()) 
        ? Channel.value([file(filtered_vcf_path), file(filtered_vcf_index_path)]) 
        : PREPARE_VCF_FILE(snps_channel, indels_channel, denylist_channel)

    // Step 10: Java installation
    log.info "Checking Java installation..."
    def java_check_result = CHECK_JAVA()
    java_check_result.java_output.view()

    // Step 11: SnpEff Tool
    log.info "Starting SnpEff setup..."

    // Check and download SnpEff tool
    def snpeff_jar_path = "${params.test_data_dir}/Tools/snpEff/snpEff.jar"
    def snpeff_config_path = "${params.test_data_dir}/Tools/snpEff/snpEff.config"
    def snpeff_tool_channel = (file(snpeff_jar_path).exists() && file(snpeff_config_path).exists()) 
        ? Channel.value([file(snpeff_jar_path), file(snpeff_config_path)]) 
        : DOWNLOAD_SNPEFF_TOOL()

    // Step 12: Check and download SnpEff database
    def snpeff_db_path = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"
    if (file(snpeff_db_path).exists()) {
        log.info "SnpEff database already exists at: ${snpeff_db_path}"
    } else {
        log.info "Downloading SnpEff database..."
        DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_tool_channel[0])
    }

    log.info "SnpEff setup completed successfully!"

    // Step 13: Arriba Tool
    def arriba_dir_path = "${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0"
    if (file(arriba_dir_path).exists()) {
        log.info "Arriba tool already exists at: ${arriba_dir_path}"
    } else {
        log.info "Downloading Arriba tool..."
        DOWNLOAD_ARRIBA()
    }

    // Step 14: VEP Cache
    def vep_cache_path = "${params.test_data_dir}/Tools/VEP/vep_cache"
    if (file(vep_cache_path).exists()) {
        log.info "VEP cache already exists at: ${vep_cache_path}"
    } else {
        log.info "Downloading VEP cache..."
        DOWNLOAD_VEP_CACHE()
    }

    // Step 15: ClinVar VCF
    def clinvar_vcf_path = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz"
    def clinvar_tbi_path = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi"
    if (file(clinvar_vcf_path).exists() && file(clinvar_tbi_path).exists()) {
        log.info "ClinVar VCF files already exist at: ${clinvar_vcf_path} and ${clinvar_tbi_path}"
    } else {
        log.info "Downloading ClinVar VCF files..."
        DOWNLOAD_CLINVAR()
    }

    // Final log
    log.info "Reference files processing workflow completed successfully!"
}

