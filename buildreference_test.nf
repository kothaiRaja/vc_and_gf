nextflow.enable.dsl = 2

// ========================== Download Reference Genome and Index ========================== //
process CHECK_OR_DOWNLOAD_REF_GENOME {
    tag "Check or Download Reference Genome"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "genome.fa", emit: genome

    when:
	!file("${params.test_data_dir}/reference/genome.fa").exists()
   
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
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

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
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

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
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

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
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output: 
    path "annotations.gtf", emit: gtf

    script:
    """
    wget -q -O annotations.gtf.gz ${params.gtf_download_url}

    # Check if the file is gzipped and unzip if necessary
    if file annotations.gtf.gz | grep -q 'gzip'; then
        gunzip annotations.gtf.gz
    else
        mv annotations.gtf.gz annotations.gtf
    fi
    """
}


// ========================== Create STAR Genome Index ========================== //
process CREATE_STAR_INDEX {
    tag "Create STAR Genome Index"
    container "https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h6b7c446_1"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

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
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed", emit: denylist

    script:
    """
    wget -q -O denylist.bed.gz ${params.denylist_download_url}

    # Check if the file is gzipped and unzip if necessary
    if file denylist.bed.gz | grep -q 'gzip'; then
        gunzip denylist.bed.gz
    else
        mv denylist.bed.gz denylist.bed
    fi
    """
}

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


// ========================== Download Indels Variants VCF ========================== //
process CHECK_OR_DOWNLOAD_VARIANTS_INDELS {
    tag "Check or Download Indels Variants"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "variants_indels.vcf.gz", emit: variants_indels

    when:
	!file("${params.test_data_dir}/reference/variants_indels.vcf.gz").exists()


    script:
    """
    wget -q -O variants_indels.vcf.gz ${params.variants_indels_download_url}
    """
}

process INDEX_INDEL_VCF {
    tag "Index INDEL VCF"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path vcf_file

    output:
    path "${vcf_file}.tbi", emit: indels_index

    script:
    """
    echo "Indexing INDEL VCF file: ${vcf_file}"
    tabix -p vcf ${vcf_file}
    """
}


process FILTER_AND_MERGE_VCF {
    tag "Filter and Merge VCF Files"
    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.15.1--h0ea216a_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

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
    publishDir "${params.test_data_dir}/Tools", mode: 'copy'
    container null

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
    tag "Download VEP Cache"
    publishDir "${params.test_data_dir}/Tools/VEP", mode: 'copy'
    container null

    output:
    path "${params.vep_cache_dir}", emit: vep_cache

    script:
    """
    mkdir -p ${params.vep_cache_dir}
    echo "Downloading VEP cache for GRCh38..."
    
    wget -q -O vep_cache.tar.gz https://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz
    tar -xzvf vep_cache.tar.gz -C ${params.vep_cache_dir} --strip-components=1
    rm vep_cache.tar.gz
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
def snpEffDbPath = ''
def known_fusions_path = ''
def blacklist_path = ''
def arribaPath = ''



workflow {

	// ========================== Define Genome Build ========================== //
	def genome_build = 'GRCh38' 
	
    // ========================== Reference Genome Handling ========================== //
    genome_ch = 
        params.reference_genome_path && file(params.reference_genome_path).exists() ? 
            Channel.of(file(params.reference_genome_path)) :
        file("${params.test_data_dir}/reference/genome.fa").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/genome.fa")) :
            CHECK_OR_DOWNLOAD_REF_GENOME()

    // Capture the reference genome path from the published directory
    genome_ch.view { genome_path ->  
    if (genome_path.toString().contains('/work/')) {
        // Force the path to the published directory
        referenceGenomePath = "${params.test_data_dir}/reference/genome.fa"
    } else {
        // If it’s already in the publish or server directory, keep it
        referenceGenomePath = genome_path.toString()
    }
    println "📂 Reference genome path set to: ${referenceGenomePath}"
}


    // ========================== Reference Genome Index Handling ========================== //
    genome_index_ch = 
        params.reference_genome_index_path && file(params.reference_genome_index_path).exists() ? 
            Channel.of(file(params.reference_genome_index_path)) :
        file("${params.test_data_dir}/reference/genome.fa.fai").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/genome.fa.fai")) :
        params.genome_index_download_url ? 
            DOWNLOAD_GENOME_INDEX() :
            CREATE_GENOME_INDEX(genome_ch)

    // Capture the genome index path from the published directory
	genome_index_ch.view { genome_index_path ->  
	if (genome_index_path.toString().contains('/work/')) {
		// Force the path to the published directory
		genomeIndexPath = "${params.test_data_dir}/reference/genome.fa.fai"
	} else {
        // If it's already in the publish or server directory, keep the original path
        genomeIndexPath = genome_index_path.toString()
    }
    println "📂 Genome index path set to: ${genomeIndexPath}"
}

	
	// ========================== Reference Genome Dictionary Handling ========================== //
    genome_dict_ch = 
        params.reference_genome_dict_path && file(params.reference_genome_dict_path).exists() ? 
            Channel.of(file(params.reference_genome_dict_path)) :
        file("${params.test_data_dir}/reference/genome.dict").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/genome.dict")) :
            CREATE_GENOME_DICT(genome_ch)

    // Capture the genome dictionary path from the published directory
	genome_dict_ch.view { genome_dict_path ->  
    if (genome_dict_path.toString().contains('/work/')) {
        // Force the path to the published directory
        genomeDictPath = "${params.test_data_dir}/reference/genome.dict"
    } else {
        // If it's already in the publish or server directory, keep the original path
        genomeDictPath = genome_dict_path.toString()
    }
    println "📂 Genome dictionary path set to: ${genomeDictPath}"
}

	
	// ========================== GTF Annotation File Handling ========================== //
    gtf_ch = 
        params.reference_genome_gtf && file(params.reference_genome_gtf).exists() ? 
            Channel.of(file(params.reference_genome_gtf)) :
        file("${params.test_data_dir}/reference/annotations.gtf").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/annotations.gtf")) :
            CHECK_OR_DOWNLOAD_GTF()

    // Capture the GTF annotation path from the published directory
gtf_ch.view { gtf_path ->  
    if (gtf_path.toString().contains('/work/')) {
        // Force the path to the published directory
        gtfPath = "${params.test_data_dir}/reference/annotations.gtf"
    } else {
        // If it's already in the publish or server directory, keep the original path
        gtfPath = gtf_path.toString()
    }
    println "📂 GTF annotation path set to: ${gtfPath}"
}

	
	// ========================== STAR Genome Index Handling ========================== //
    star_index_ch = 
        params.star_genome_index_path && file(params.star_genome_index_path).exists() ? 
            Channel.of(file(params.star_genome_index_path)) :
        file("${params.test_data_dir}/reference/STAR_index").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/STAR_index")) :
            CREATE_STAR_INDEX(genome_ch, gtf_ch)

    // Capture the STAR genome index path from the published directory
star_index_ch.view { star_index_path ->  
    if (star_index_path.toString().contains('/work/')) {
        // Force the path to the published directory
        starIndexPath = "${params.test_data_dir}/reference/STAR_index"
    } else {
        // If it's already in the publish or server directory, keep the original path
        starIndexPath = star_index_path.toString()
    }
    println "📂 STAR genome index path set to: ${starIndexPath}"
}

	
	// ========================== Denylist File Handling ========================== //
    denylist_ch = 
        params.reference_denylist_path && file(params.reference_denylist_path).exists() ? 
            Channel.of(file(params.reference_denylist_path)) :
        file("${params.test_data_dir}/reference/denylist.bed").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/denylist.bed")) :
            CHECK_OR_DOWNLOAD_DENYLIST()

    // Capture the denylist path from the published directory
denylist_ch.view { denylist_path ->  
    if (denylist_path.toString().contains('/work/')) {
        // Force the path to the published directory
        denylistPath = "${params.test_data_dir}/reference/denylist.bed"
    } else {
        // If it's already in the publish or server directory, keep the original path
        denylistPath = denylist_path.toString()
    }
    println "📂 Denylist BED file path set to: ${denylistPath}"
}

	
	// ========================== SNP VCF Handling ========================== //
    snp_vcf_ch = 
        params.variants_snp_path && file(params.variants_snp_path).exists() ? 
            Channel.of(file(params.variants_snp_path)) :
        file("${params.test_data_dir}/reference/variants_snp.vcf.gz").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/variants_snp.vcf.gz")) :
            CHECK_OR_DOWNLOAD_VARIANTS_SNP()
	
	//Capture the snps path from published directory
	
	snp_vcf_ch.view { snp_vcf_path ->  
    if (snp_vcf_path.toString().contains('/work/')) {
        // Force the path to the published directory
        snpVcfPath = "${params.test_data_dir}/reference/variants_snp.vcf.gz"
    } else {
        // If it's already in the publish or server directory, keep the original path
        snpVcfPath = snp_vcf_path.toString()
    }
    println "📂 SNP VCF path set to: ${snpVcfPath}"
}


    // ========================== SNP Index Handling ========================== //
    snp_index_ch = 
    params.variants_snp_index_path && file(params.variants_snp_index_path).exists() ? 
        Channel.of(file(params.variants_snp_index_path)) :
    file("${params.test_data_dir}/reference/variants_snp.vcf.gz.tbi").exists() ? 
        Channel.of(file("${params.test_data_dir}/reference/variants_snp.vcf.gz.tbi")) :
    params.variants_snp_index_download_url ? 
        DOWNLOAD_VARIANTS_SNP_INDEX() :
        INDEX_SNP_VCF(snp_vcf_ch)
			
	snp_index_ch.view { snp_index_path ->  
    if (snp_index_path.toString().contains('/work/')) {
        // Force the path to the published directory
        snpIndexPath = "${params.test_data_dir}/reference/variants_snp.vcf.gz.tbi"
    } else {
        // If it's already in the publish or server directory, keep the original path
        snpIndexPath = snp_index_path.toString()
    }
    println "📂 SNP Index path set to: ${snpIndexPath}"
}


	// ========================== Indels VCF Handling ========================== //

    indels_vcf_ch = 
        params.variants_indels_path && file(params.variants_indels_path).exists() ? 
            Channel.of(file(params.variants_indels_path)) :
        file("${params.test_data_dir}/reference/variants_indels.vcf.gz").exists() ?
            Channel.of(file("${params.test_data_dir}/reference/variants_indels.vcf.gz")) :
            CHECK_OR_DOWNLOAD_VARIANTS_INDELS()
			
	indels_vcf_ch.view { indels_vcf_path ->  
    if (indels_vcf_path.toString().contains('/work/')) {
        // Force the path to the published directory
        indelsVcfPath = "${params.test_data_dir}/reference/variants_indels.vcf.gz"
    } else {
        // If it's already in the publish or server directory, keep the original path
        indelsVcfPath = indels_vcf_path.toString()
    }
    println "📂 Indels VCF path set to: ${indelsVcfPath}"
}

	// ========================== Indels Index Handling ========================== //

	indels_index_ch = 
    params.variants_indels_index_path && file(params.variants_indels_index_path).exists() ? 
        Channel.of(file(params.variants_indels_index_path)) :
    file("${params.test_data_dir}/reference/variants_indels.vcf.gz.tbi").exists() ? 
        Channel.of(file("${params.test_data_dir}/reference/variants_indels.vcf.gz.tbi")) :
    params.variants_indels_index_download_url ? 
        DOWNLOAD_VARIANTS_INDELS_INDEX() :
        INDEX_INDEL_VCF(indels_vcf_ch).indels_index
		
	indels_index_ch.view { indels_index_path ->  
    if (indels_index_path.toString().contains('/work/')) {
        // Redirect to the published directory
        indelsIndexPath = "${params.test_data_dir}/reference/variants_indels.vcf.gz.tbi"
    } else {
        // Keep the original path if it's from the server or already published
        indelsIndexPath = indels_index_path.toString()
    }
    println "📂 Indels Index path set to: ${indelsIndexPath}"
}


	// ========================== Filter and Merge VCFs ========================== //
    
	
// Define channels for merged VCF and index
    def merged_vcf_ch
    def merged_vcf_index_ch

    if (file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz").exists() &&
        file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi").exists()) {

        // If files exist in the publish directory, create channels from them
        merged_vcf_ch = Channel.of(file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz"))
        merged_vcf_index_ch = Channel.of(file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"))

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
    if (merged_vcf_path.toString().contains('/work/')) {
        // Redirect to the published directory
        mergedVcfPath = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz"
    } else {
        // Keep the original path if it's from the server or already published
        mergedVcfPath = merged_vcf_path.toString()
    }
    println "📂 Merged VCF path set to: ${mergedVcfPath}"
}

// Capture merged VCF index path
merged_vcf_index_ch.view { merged_vcf_index_path ->
    if (merged_vcf_index_path.toString().contains('/work/')) {
        // Redirect to the published directory
        mergedVcfIndexPath = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"
    } else {
        // Keep the original path if it's from the server or already published
        mergedVcfIndexPath = merged_vcf_index_path.toString()
    }
    println "📂 Merged VCF Index path set to: ${mergedVcfIndexPath}"
}

	
	// ========================== Java Version Check ========================== //

    def java_check_ch

    if (file("${params.test_data_dir}/reference/java_check.log").exists()) {
        println "✅ Java check log already exists in the publish directory. Skipping Java check."

        // If the log exists, create a channel from the existing log file
        java_check_ch = Channel.of(file("${params.test_data_dir}/reference/java_check.log"))

    } else {
        // Run the CHECK_JAVA process if the log file does not exist
        java_check_ch = CHECK_JAVA().java_output
    }

    // ========================== Capture Java Check Output ========================== //

    java_check_ch.view { java_log_path ->
        javaLogPath = file("${params.test_data_dir}/reference/java_check.log").toString()
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

	} else if (file("${params.test_data_dir}/Tools/snpEff.jar").exists() && 
			file("${params.test_data_dir}/Tools/snpEff.config").exists()) {
    
		println "✅ SnpEff tool found in the publish directory. Skipping download."

		snpeff_jar_ch = Channel.of(file("${params.test_data_dir}/Tools/snpEff.jar"))
		snpeff_config_ch = Channel.of(file("${params.test_data_dir}/Tools/snpEff.config"))

		snpEffJarPath = "${params.test_data_dir}/Tools/snpEff.jar"
		snpEffConfigPath = "${params.test_data_dir}/Tools/snpEff.config"

	} else {
		println "⚠️ SnpEff tool not found. Downloading..."
		def result = DOWNLOAD_SNPEFF_TOOL()
		snpeff_jar_ch = result.snpeff_jar
		snpeff_config_ch = result.snpeff_config
}


    // ========================== Capture SnpEff Paths ========================== //

    // Capture SnpEff JAR path
snpeff_jar_ch.view { snpeff_jar_path ->  
    if (snpeff_jar_path.toString().contains('/work/')) {
        // Redirect to the published directory
        snpEffJarPath = "${params.test_data_dir}/Tools/snpEff/snpEff.jar"
    } else {
        // Keep the original path if it's from the server or already published
        snpEffJarPath = snpeff_jar_path.toString()
    }
    println "📂 SnpEff JAR path set to: ${snpEffJarPath}"
}

// Capture SnpEff Config path
snpeff_config_ch.view { snpeff_config_path ->  
    if (snpeff_config_path.toString().contains('/work/')) {
        // Redirect to the published directory
        snpEffConfigPath = "${params.test_data_dir}/Tools/snpEff/snpEff.config"
    } else {
        // Keep the original path if it's from the server or already published
        snpEffConfigPath = snpeff_config_path.toString()
    }
    println "📂 SnpEff Config path set to: ${snpEffConfigPath}"
}

	
	// ========================== SnpEff Database Handling ========================== //

	def snpeff_db_ch  // Channel to handle the database path dynamically

if (params.snpeff_db_dir && file("${params.snpeff_db_dir_path}/${params.genomedb}").exists()) {
    println "✅ SnpEff database for ${params.genomedb} found in the server directory. Skipping download."

    snpeff_db_ch = Channel.of(file("${params.snpeff_db_dir_path}/${params.genomedb}"))
    snpEffDbPath = "${params.snpeff_db_dir_path}/${params.genomedb}"

} else if (file("${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}").exists()) {
    println "✅ SnpEff database found in the publish directory. Skipping download."

    snpeff_db_ch = Channel.of(file("${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"))
    snpEffDbPath = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"

} else {
    println "⚠️ SnpEff database not found. Downloading..."
    def result = DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_jar_ch)

    snpeff_db_ch = result
    snpEffDbPath = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"
}

// ========================== Capture SnpEff Database Path ========================== //

// Capture SnpEff Database path
snpeff_db_ch.view { snpeff_db_path ->  
    if (snpeff_db_path.toString().contains('/work/')) {
        // Redirect to the published directory
        snpEffDbPath = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"
    } else {
        // Keep the original path if it's from the server or already published
        snpEffDbPath = snpeff_db_path.toString()
    }
    println "📂 SnpEff Database path set to: ${snpEffDbPath}"
}



// ========================== Arriba Tool Handling ========================== //


def arriba_dir_ch

if (params.arriba_tool_dir_path && file("${params.arriba_tool_dir_path}/arriba_v2.4.0").exists()) {
    println "✅ Arriba tool found in the server directory."
    arriba_dir_ch = Channel.of(file("${params.arriba_tool_dir_path}/arriba_v2.4.0"))

} else if (file("${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0").exists()) {
    println "✅ Arriba tool found in the publish directory."
    arriba_dir_ch = Channel.of(file("${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0"))

} else {
    println "⚠️ Arriba tool not found. Downloading..."
    def result = DOWNLOAD_ARRIBA()
    arriba_dir_ch = result.arriba_dir
}

// ========================== Capture Arriba Tool and Database Paths ========================== //
// ========================== Handle Arriba Tool and Capture Paths ========================== //
arriba_dir_ch.view { arriba_dir_path ->  
    // Set Arriba path after checking if it exists in the server or publish directory
    if (arriba_dir_path.toString().contains('/work/')) {
        // If downloaded (i.e., in workdir), set to publish directory
        arribaPath = "${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0"
    } else {
        // If from server or already published, use the existing path
        arribaPath = arriba_dir_path.toString()
    }

    println "📂 Arriba tool path set to: ${arribaPath}"

    // Define known fusions and blacklist paths directly without immediate file check
    knownFusionsPath = "${arribaPath}/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
    blacklistPath = "${arribaPath}/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"

    println "📂 Known fusions path set to: ${knownFusionsPath}"
    println "📂 Blacklist path set to: ${blacklistPath}"
}




// ========================== VEP Cache Handling ========================== //
    def vep_cache_ch

    if (params.vep_cache_dir_path && file("${params.vep_cache_dir_path}").exists()) {
        println "✅ VEP cache found in the server directory."
        vep_cache_ch = Channel.of(file("${params.vep_cache_dir_path}"))

    } else if (file("${params.test_data_dir}/Tools/VEP").exists()) {
        println "✅ VEP cache found in the publish directory."
        vep_cache_ch = Channel.of(file("${params.test_data_dir}/Tools/VEP"))

    } else {
        println "⚠️ VEP cache not found. Downloading..."
        def result = DOWNLOAD_VEP_CACHE()
        vep_cache_ch = result.vep_cache
    }

    // ========================== Capture VEP Cache Path ========================== //
vep_cache_ch.view { vep_cache_path ->  
    if (vep_cache_path.toString().contains('/work/')) {
        // Redirect to the published directory
        vepCachePath = "${params.test_data_dir}/Tools/VEP"
    } else {
        // Keep the original path if it's from the server or already published
        vepCachePath = vep_cache_path.toString()
    }
    println "📂 VEP Cache path set to: ${vepCachePath}"
}


// ========================== ClinVar VCF Handling ========================== //
    def clinvar_vcf_ch, clinvar_tbi_ch

    if (file("${params.clinvar_path}").exists() && file("${params.clinvartbi_path}").exists()) {
        println "✅ ClinVar VCF and index found in the server directory."
        clinvar_vcf_ch = Channel.of(file("${params.clinvar_path}"))
        clinvar_tbi_ch = Channel.of(file("${params.clinvartbi_path}"))

    } else if (file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz").exists() && 
               file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi").exists()) {
        println "✅ ClinVar VCF and index found in the publish directory."
        clinvar_vcf_ch = Channel.of(file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz"))
        clinvar_tbi_ch = Channel.of(file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi"))

    } else {
        println "⚠️ ClinVar VCF and index not found. Downloading..."
        def result = DOWNLOAD_CLINVAR()
        clinvar_vcf_ch = result.out[0]
        clinvar_tbi_ch = result.out[1]
    }

    // ========================== Capture ClinVar Paths ========================== //
    clinvar_vcf_ch.view { clinvar_vcf_path ->  
    if (clinvar_vcf_path.toString().contains('/work/')) {
        // Redirect to the published directory
        clinvarVcfPath = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz"
    } else {
        // Keep the original path if it's from the server or already published
        clinvarVcfPath = clinvar_vcf_path.toString()
    }
    println "📂 ClinVar VCF path set to: ${clinvarVcfPath}"
}

clinvar_tbi_ch.view { clinvar_tbi_path ->  
    if (clinvar_tbi_path.toString().contains('/work/')) {
        // Redirect to the published directory
        clinvarTbiPath = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi"
    } else {
        // Keep the original path if it's from the server or already published
        clinvarTbiPath = clinvar_tbi_path.toString()
    }
    println "📂 ClinVar VCF Index path set to: ${clinvarTbiPath}"
}





    // ========================== Writing to Config After Completion ========================== //
    workflow.onComplete {
        try {
            def baseDir = System.getProperty('user.dir')
            def outputDir = new File("${baseDir}/reference_test_paths.config").getParentFile()

            if (!outputDir.exists()) {
                println "📁 Output directory does not exist. Creating: ${outputDir}"
                outputDir.mkdirs()
            }

            def configFile = new File("${baseDir}/reference_test_paths.config")

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
			params.arriba_tool_dir = '${arribaPath ?: 'NOT_FOUND'}'
			params.arriba_known_fusions = '${knownFusionsPath ?: 'NOT_FOUND'}'
			params.arriba_blacklist = '${blacklistPath ?: 'NOT_FOUND'}'
			params.vep_cache_dir = '${vepCachePath ?: 'NOT_FOUND'}'
			params.clinvar = '${clinvarVcfPath ?: 'NOT_FOUND'}'
            params.clinvartbi = '${clinvarTbiPath ?: 'NOT_FOUND'}'
            """

            println "✅ Reference paths successfully written to ${configFile}"

        } catch (Exception e) {
            println "❌ Error writing reference paths: ${e.message}"
        }
    }
}
