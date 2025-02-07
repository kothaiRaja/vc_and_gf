nextflow.enable.dsl = 2

// ========================== Download Reference Genome and Index ========================== //
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

process DOWNLOAD_GENOME_INDEX {
    tag "Download Genome Index"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa.fai", emit: genome_fai

    script:
    """
    echo "⚠️ Downloading genome index from provided URL..."
    wget -q -O genome.fa.fai ${params.genome_index_download_url}
    """
}

process CREATE_GENOME_INDEX {
    tag "Create Genome Index"
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.14--hb421002_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa

    output:
    path "genome.fa.fai", emit: genome_fai

    script:
    """
    echo "⚠️ Creating genome index using samtools..."
    samtools faidx $genome_fa
    """
}

// ========================== Create Genome Dictionary ========================== //
process CREATE_GENOME_DICT {
    tag "Create Genome Dictionary"
    container "https://depot.galaxyproject.org/singularity/picard%3A2.27.4--hdfd78af_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa

    output:
    path "genome.dict", emit: genome_dict

    script:
    """
    echo "⚠️ Creating genome dictionary using Picard..."
    picard CreateSequenceDictionary R=$genome_fa O=genome.dict
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


// ========================== Create STAR Genome Index ========================== //
process CREATE_STAR_INDEX {
    tag "Create STAR Genome Index"
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

// ========================== Download SNP Variants Index ========================== //
process DOWNLOAD_VARIANTS_SNP_INDEX {
    tag "Download SNP Index"
    container null
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    output:
    path "variants_snp.vcf.gz.tbi", emit: snp_index

    script:
    """
    wget -q -O variants_snp.vcf.gz.tbi ${params.variants_snp_index_download_url}
    """
}

process INDEX_VCF {
    tag "Index VCF File"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input:
    path vcf_file

    output:
    path "${vcf_file}.tbi", emit: vcf_index

    script:
    """
    echo \"⚠️ Indexing VCF file using tabix...\"
    tabix -p vcf ${vcf_file}
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

process FILTER_AND_MERGE_VCF {
    tag "Filter and Merge VCF Files"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.actual_data_dir}/reference", mode: 'copy'

    input: 
    path variants_snp
    path variants_snp_index 
    path variants_indels
    path variants_indels_index   
    path denylist

    output:
    path "merged.filtered.recode.vcf.gz", emit: merged_vcf
    path "merged.filtered.recode.vcf.gz.tbi", emit: merged_vcf_tbi

    script:
    """
    # Filter SNPs (bcftools will automatically use .tbi if it exists)
    bcftools view -T ^${denylist} ${variants_snp} -Oz -o filtered_snps.vcf.gz
    tabix -p vcf filtered_snps.vcf.gz  # Ensure indexing if not present

    # Filter INDELs
    bcftools view -T ^${denylist} ${variants_indels} -Oz -o filtered_indels.vcf.gz
    tabix -p vcf filtered_indels.vcf.gz  # Ensure indexing if not present

    # Merge the filtered SNPs and INDELs into a single VCF file
    bcftools merge filtered_snps.vcf.gz filtered_indels.vcf.gz -Oz -o merged.filtered.recode.vcf.gz
    tabix -p vcf merged.filtered.recode.vcf.gz  # Index the merged VCF

    # Check and fix contig prefixes before finalizing
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
    container null  // No container needed for local Java check
	
    output:
    path "java_check.log", emit: java_output
	
    when:
    true  // Always run to check Java status

    script:
    """
    # Get the Java version
    java_version=\$(java -version 2>&1 | grep -oP '(?<=version ")([0-9]+)' | head -1)
    
    if [[ -z "\$java_version" ]]; then
        echo "❌ ERROR: Java is not installed or not in PATH."
        echo "❌ ERROR: Java is not installed or not in PATH." > java_check.log
    elif [[ \$java_version -lt 21 ]]; then
        echo "⚠️ WARNING: Java version \$java_version detected. Please update to version 21 or higher."
        echo "⚠️ WARNING: Java version \$java_version detected. Please update to version 21 or higher." > java_check.log
    else
        echo "✅ Java version \$java_version detected. Java is up-to-date."
        echo "✅ Java version \$java_version detected. Java is up-to-date." > java_check.log
    fi
    """
}

process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
    publishDir "${params.actual_data_dir}/Tools", mode: 'copy'
    container null  // No container needed, using wget and unzip

    output:
    path "${params.snpeff_jar_dir}/snpEff.jar", emit: snpeff_jar
    path "${params.snpeff_jar_dir}/snpEff.config", emit: snpeff_config

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


// Define variables to store the reference genome and genome index paths
def referenceGenomePath = ''
def genomeIndexPath = ''
def genomeDictPath = ''
def gtfPath = ''
def starIndexPath = ''
def snpVcfPath = ''
def snpIndexPath = '' 
def indelsVcfPath = ''
def indelsIndexPath = ''
def mergedVcfPath = ''
def mergedVcfIndexPath = ''



workflow {
    // ========================== Reference Genome Handling ========================== //
    genome_ch = 
        params.reference_genome_path && file(params.reference_genome_path).exists() ? 
            Channel.of(file(params.reference_genome_path)) :
        file("${params.actual_data_dir}/reference/genome.fa").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/genome.fa")) :
            CHECK_OR_DOWNLOAD_REF_GENOME()

    // Capture the reference genome path from the published directory
    genome_ch.view { path ->
        referenceGenomePath = file("${params.actual_data_dir}/reference/genome.fa").toString()
        println "📂 Reference genome path set to: ${referenceGenomePath}"
    }

    // ========================== Reference Genome Index Handling ========================== //
    genome_index_ch = 
        params.reference_genome_index_path && file(params.reference_genome_index_path).exists() ? 
            Channel.of(file(params.reference_genome_index_path)) :
        file("${params.actual_data_dir}/reference/genome.fa.fai").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/genome.fa.fai")) :
        params.genome_index_download_url ? 
            DOWNLOAD_GENOME_INDEX.out.genome_fai :
            CREATE_GENOME_INDEX(genome_ch)

    // Capture the genome index path from the published directory
    genome_index_ch.view { path ->
        genomeIndexPath = file("${params.actual_data_dir}/reference/genome.fa.fai").toString()
        println "📂 Genome index path set to: ${genomeIndexPath}"
    }
	
	// ========================== Reference Genome Dictionary Handling ========================== //
    genome_dict_ch = 
        params.reference_genome_dict_path && file(params.reference_genome_dict_path).exists() ? 
            Channel.of(file(params.reference_genome_dict_path)) :
        file("${params.actual_data_dir}/reference/genome.dict").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/genome.dict")) :
            CREATE_GENOME_DICT(genome_ch)

    // Capture the genome dictionary path from the published directory
    genome_dict_ch.view { path ->
        genomeDictPath = file("${params.actual_data_dir}/reference/genome.dict").toString()
        println "📂 Genome dictionary path set to: ${genomeDictPath}"
    }
	
	// ========================== GTF Annotation File Handling ========================== //
    gtf_ch = 
        params.reference_genome_gtf && file(params.reference_genome_gtf).exists() ? 
            Channel.of(file(params.reference_genome_gtf)) :
        file("${params.actual_data_dir}/reference/annotations.gtf").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/annotations.gtf")) :
            CHECK_OR_DOWNLOAD_GTF()

    // Capture the GTF annotation path from the published directory
    gtf_ch.view { path ->
        gtfPath = file("${params.actual_data_dir}/reference/annotations.gtf").toString()
        println "📂 GTF annotation path set to: ${gtfPath}"
    }
	
	// ========================== STAR Genome Index Handling ========================== //
    star_index_ch = 
        params.star_genome_index_path && file(params.star_genome_index_path).exists() ? 
            Channel.of(file(params.star_genome_index_path)) :
        file("${params.actual_data_dir}/reference/STAR_index").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/STAR_index")) :
            CREATE_STAR_INDEX(genome_ch, gtf_ch)

    // Capture the STAR genome index path from the published directory
    star_index_ch.view { path ->
        starIndexPath = file("${params.actual_data_dir}/reference/STAR_index").toString()
        println "📂 STAR genome index path set to: ${starIndexPath}"
    }
	
	// ========================== Denylist File Handling ========================== //
    denylist_ch = 
        params.reference_denylist_path && file(params.reference_denylist_path).exists() ? 
            Channel.of(file(params.reference_denylist_path)) :
        file("${params.actual_data_dir}/reference/denylist.bed").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/denylist.bed")) :
            CHECK_OR_DOWNLOAD_DENYLIST()

    // Capture the denylist path from the published directory
    denylist_ch.view { path ->
        denylistPath = file("${params.actual_data_dir}/reference/denylist.bed").toString()
        println "📂 Denylist BED file path set to: ${denylistPath}"
    }
	
	// ========================== SNP VCF Handling ========================== //
    snp_vcf_ch = 
        params.variants_snp_path && file(params.variants_snp_path).exists() ? 
            Channel.of(file(params.variants_snp_path)) :
        file("${params.actual_data_dir}/reference/variants_snp.vcf.gz").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/variants_snp.vcf.gz")) :
            CHECK_OR_DOWNLOAD_VARIANTS_SNP()
	
	//Capture the snps path from published directory
	
	snp_vcf_ch.view { snp_vcf_path ->  // Changed variable name to 'snp_vcf_path'
    snpVcfPath = file("${params.actual_data_dir}/reference/variants_snp.vcf.gz").toString()
    println "📂 SNP VCF path set to: ${snpVcfPath}"
}

    // ========================== SNP Index Handling ========================== //
    snp_index_ch = 
        params.variants_snp_index_path && file(params.variants_snp_index_path).exists() ? 
            Channel.of(file(params.variants_snp_index_path)) :
        file("${params.actual_data_dir}/reference/variants_snp.vcf.gz.tbi").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/variants_snp.vcf.gz.tbi")) :
        params.variants_snp_index_download_url ? 
            DOWNLOAD_VARIANTS_SNP_INDEX() :
            INDEX_VCF(snp_vcf_ch)
			
	// Capture the SNP Index path from the published directory
	snp_index_ch.view { snp_index_path ->  
    snpIndexPath = file("${params.actual_data_dir}/reference/variants_snp.vcf.gz.tbi").toString()
    println "📂 SNP Index path set to: ${snpIndexPath}"
}

	// ========================== Indels VCF Handling ========================== //

    indels_vcf_ch = 
        params.variants_indels_path && file(params.variants_indels_path).exists() ? 
            Channel.of(file(params.variants_indels_path)) :
        file("${params.actual_data_dir}/reference/variants_indels.vcf.gz").exists() ?
            Channel.of(file("${params.actual_data_dir}/reference/variants_indels.vcf.gz")) :
            CHECK_OR_DOWNLOAD_VARIANTS_INDELS()
			
	//Capture the indels path from published directory
	
	indels_vcf_ch.view { indels_vcf_path ->  
    indelsVcfPath = file("${params.actual_data_dir}/reference/variants_indels.vcf.gz").toString()
    println "📂 Indels VCF path set to: ${indelsVcfPath}"
}

	// ========================== Indels Index Handling ========================== //

	indels_index_ch = 
    params.variants_indels_index_path && file(params.variants_indels_index_path).exists() ? 
        Channel.of(file(params.variants_indels_index_path)) :
    file("${params.actual_data_dir}/reference/variants_indels.vcf.gz.tbi").exists() ?
        Channel.of(file("${params.actual_data_dir}/reference/variants_indels.vcf.gz.tbi")) :
    params.variants_indels_index_download_url ? 
        DOWNLOAD_VARIANTS_INDELS_INDEX() :
        INDEX_VCF(indels_vcf_ch) 
		
	// Capture the Indels Index path from the published directory
	indels_index_ch.view { indels_index_path ->  
    indelsIndexPath = file("${params.actual_data_dir}/reference/variants_indels.vcf.gz.tbi").toString()
    println "📂 Indels Index path set to: ${indelsIndexPath}"
}

	// ========================== Filter and Merge VCFs ========================== //
    
	
// Define channels for merged VCF and index
    def merged_vcf_ch
    def merged_vcf_index_ch

    if (file("${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz").exists() &&
        file("${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi").exists()) {

        // If files exist in the publish directory, create channels from them
        merged_vcf_ch = Channel.of(file("${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz"))
        merged_vcf_index_ch = Channel.of(file("${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"))

        println "✅ Merged VCF and index already exist in the publish directory. Skipping merge."

    } else {
        // Run the FILTER_AND_MERGE_VCF process if files don't exist
         merge_result = FILTER_AND_MERGE_VCF(snp_vcf_ch, snp_index_ch, indels_vcf_ch, indels_index_ch, denylist_ch)

        // Capture the outputs from the process
        merged_vcf_ch = merge_result.merged_vcf
        merged_vcf_index_ch = merge_result.merged_vcf_tbi
    }

    // ========================== Capture Merged VCF Paths ========================== //

    // Capture merged VCF path
    merged_vcf_ch.view { merged_vcf_path ->
        mergedVcfPath = file("${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz").toString()
        println "📂 Merged VCF path set to: ${mergedVcfPath}"
    }

    // Capture merged VCF index path
    merged_vcf_index_ch.view { merged_vcf_index_path ->
        mergedVcfIndexPath = file("${params.actual_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi").toString()
        println "📂 Merged VCF Index path set to: ${mergedVcfIndexPath}"
    }
	
	// ========================== Java Version Check ========================== //

    def java_check_ch

    if (file("${params.actual_data_dir}/reference/java_check.log").exists()) {
        println "✅ Java check log already exists in the publish directory. Skipping Java check."

        // If the log exists, create a channel from the existing log file
        java_check_ch = Channel.of(file("${params.actual_data_dir}/reference/java_check.log"))

    } else {
        // Run the CHECK_JAVA process if the log file does not exist
        java_check_ch = CHECK_JAVA().java_output
    }

    // ========================== Capture Java Check Output ========================== //

    java_check_ch.view { java_log_path ->
        javaLogPath = file("${params.actual_data_dir}/reference/java_check.log").toString()
        println "📂 Java check log path set to: ${javaLogPath}"
    }
	
	// ========================== SnpEff Tool Handling ========================== //

    def snpeff_jar_ch, snpeff_config_ch

	if (file(params.snpeff_jar_path).exists() && file(params.snpeff_config_path).exists()) {
		println "✅ SnpEff tool found in the server directory. Skipping download."

		snpeff_jar_ch = Channel.of(file(params.snpeff_jar_path))
		snpeff_config_ch = Channel.of(file(params.snpeff_config_path))

		snpEffJarPath = params.snpeff_jar_path
		snpEffConfigPath = params.snpeff_config_path

	} else if (file("${params.actual_data_dir}/Tools/snpEff.jar").exists() && 
			file("${params.actual_data_dir}/Tools/snpEff.config").exists()) {
    
		println "✅ SnpEff tool found in the publish directory. Skipping download."

		snpeff_jar_ch = Channel.of(file("${params.actual_data_dir}/Tools/snpEff.jar"))
		snpeff_config_ch = Channel.of(file("${params.actual_data_dir}/Tools/snpEff.config"))

		snpEffJarPath = "${params.actual_data_dir}/Tools/snpEff.jar"
		snpEffConfigPath = "${params.actual_data_dir}/Tools/snpEff.config"

	} else {
		println "⚠️ SnpEff tool not found. Downloading..."
		def result = DOWNLOAD_SNPEFF_TOOL()
		snpeff_jar_ch = result.snpeff_jar
		snpeff_config_ch = result.snpeff_config
}


    // ========================== Capture SnpEff Paths ========================== //

    snpeff_jar_ch.view { snpeff_jar_path ->  
        snpEffJarPath = ("${params.actual_data_dir}/Tools/snpEff.jar").toString()
        println "📂 SnpEff JAR path set to: ${snpEffJarPath}"
    }

    snpeff_config_ch.view { snpeff_config_path ->  
        snpEffConfigPath = ("${params.actual_data_dir}/Tools/snpEff.config").toString()
        println "📂 SnpEff Config path set to: ${snpEffConfigPath}"
    }
	
	// ========================== SnpEff Database Handling ========================== //

	def snpeff_db_ch  // Channel to handle the database path dynamically

	if (params.snpeff_db_dir && file("${params.snpeff_db_dir_path}/${params.genomedb}").exists()) {
		println "✅ SnpEff database for ${params.genomedb} found in the server directory. Skipping download."

		snpeff_db_ch = Channel.of(file("${params.snpeff_db_dir_path}/${params.genomedb}"))
		snpEffDbPath = "${params.snpeff_db_dir}/${params.genomedb}"

	} else if (file("${params.actual_data_dir}/Tools/snpEff/${params.genomedb}").exists()) {
		println "✅ SnpEff database found in the publish directory. Skipping download."

		snpeff_db_ch = Channel.of(file("${params.actual_data_dir}/Tools/snpEff/${params.genomedb}"))
		snpEffDbPath = "${params.actual_data_dir}/Tools/snpEff/${params.genomedb}"

	} else {
		println "⚠️ SnpEff database not found. Downloading..."
		def result = DOWNLOAD_SNPEFF_DB(genome: params.genomedb, snpeff_jar_path: snpeff_jar_ch)

		snpeff_db_ch = result
		snpEffDbPath = "${params.actual_data_dir}/Tools/snpEff/${params.genomedb}"
}

// ========================== Capture SnpEff Database Path ========================== //

	snpeff_db_ch.view { snpeff_db_path ->  
		snpEffDbPath = "${params.actual_data_dir}/Tools/snpEff/${params.genomedb}".toString()
		println "📂 SnpEff Database path set to: ${snpEffDbPath}"
}

// ========================== Arriba Tool Handling ========================== //

def arriba_dir_ch  // Channel to dynamically handle Arriba's path

if (params.arriba_tool_dir && file("${params.arriba_tool_dir}/arriba_v2.4.0").exists()) {
    println "✅ Arriba tool found in the server directory. Skipping download."

    arriba_dir_ch = Channel.of(file("${params.arriba_tool_dir}/arriba_v2.4.0"))
    arribaPath = "${params.arriba_tool_dir}/arriba_v2.4.0"

} else if (file("${params.actual_data_dir}/Tools/ARRIBA/arriba_v2.4.0").exists()) {
    println "✅ Arriba tool found in the publish directory. Skipping download."

    arriba_dir_ch = Channel.of(file("${params.actual_data_dir}/Tools/ARRIBA/arriba_v2.4.0"))
    arribaPath = "${params.actual_data_dir}/Tools/ARRIBA/arriba_v2.4.0"

} else {
    println "⚠️ Arriba tool not found. Downloading..."
    def result = DOWNLOAD_ARRIBA()

    arriba_dir_ch = result.arriba_dir
    arribaPath = "${params.actual_data_dir}/Tools/ARRIBA/arriba_v2.4.0"
}

// ========================== Capture Arriba Path ========================== //

arriba_dir_ch.view { path ->  
	arribaPath = "${params.actual_data_dir}/Tools/ARRIBA".toString()
    println "📂 Arriba tool path set to: ${arribaPath}"
}





    // ========================== Writing to Config After Completion ========================== //
    workflow.onComplete {
        try {
            def baseDir = System.getProperty('user.dir')
            def outputDir = new File("${baseDir}/reference_paths.config").getParentFile()

            if (!outputDir.exists()) {
                println "📁 Output directory does not exist. Creating: ${outputDir}"
                outputDir.mkdirs()
            }

            def configFile = new File("${baseDir}/reference_paths.config")

            configFile.text = """  
            params.reference_genome = '${referenceGenomePath ?: 'NOT_FOUND'}'
            params.reference_genome_index = '${genomeIndexPath ?: 'NOT_FOUND'}'
			params.reference_genome_dict = '${genomeDictPath ?: 'NOT_FOUND'}'
			params.gtf_annotation = '${gtfPath ?: 'NOT_FOUND'}'
			params.star_genome_index = '${starIndexPath ?: 'NOT_FOUND'}'
			params.denylist_bed = '${denylistPath ?: 'NOT_FOUND'}'
			params.variants_snp = '${snpVcfPath ?: 'NOT_FOUND'}'
            params.variants_snp_index = '${snpIndexPath ?: 'NOT_FOUND'}'
			params.variants_indels = '${indelsVcfPath ?: 'NOT_FOUND'}'
            params.variants_indels_index = '${indelsIndexPath ?: 'NOT_FOUND'}'
			params.merged_vcf = '${mergedVcfPath ?: 'NOT_FOUND'}'
			params.merged_vcf_index = '${mergedVcfIndexPath ?: 'NOT_FOUND'}'
			params.snpeff_jar = '${snpEffJarPath ?: 'NOT_FOUND'}'
            params.snpeff_config = '${snpEffConfigPath ?: 'NOT_FOUND'}'
			params.snpeff_db = '${snpEffDbPath ?: 'NOT_FOUND'}'
			params.arriba_tool = '${arribaPath ?: 'NOT_FOUND'}'
			
            """

            println "✅ Reference paths successfully written to ${configFile}"

        } catch (Exception e) {
            println "❌ Error writing reference paths: ${e.message}"
        }
    }
}
