nextflow.enable.dsl = 2

include { CHECK_OR_DOWNLOAD_REF_GENOME } from './modules/references/Test/reference_genome_test.nf'
include { DOWNLOAD_GENOME_INDEX } from './modules/references/Test/genome_index_test.nf'
include { CREATE_GENOME_INDEX } from './modules/references/Test/genome_index_test.nf'
include { CREATE_GENOME_DICT } from './modules/references/Test/genome_dict_test.nf'
include { CHECK_OR_DOWNLOAD_GTF } from './modules/references/Test/gtf_annotation_test.nf'
include { CREATE_STAR_INDEX } from './modules/references/Test/STAR_index_test.nf'
include { CHECK_OR_DOWNLOAD_DENYLIST } from './modules/references/Test/denylist_test.nf'
include { CHECK_OR_DOWNLOAD_VARIANTS_SNP } from './modules/references/Test/snp_test.nf'
include { DOWNLOAD_VARIANTS_SNP_INDEX } from './modules/references/Test/snp_test.nf'
include { INDEX_SNP_VCF } from './modules/references/Test/snp_test.nf'
include { CHECK_OR_DOWNLOAD_VARIANTS_INDELS } from './modules/references/Test/indels_test.nf'
include { INDEX_INDEL_VCF } from './modules/references/Test/indels_test.nf'
include { DOWNLOAD_VARIANTS_INDELS_INDEX } from './modules/references/Test/indels_test.nf'
include { FILTER_AND_MERGE_VCF } from './modules/references/Test/prepare_vcfs_test.nf'
include { CHECK_JAVA } from './modules/references/check_java.nf'
include { DOWNLOAD_SNPEFF_TOOL } from './modules/references/Test/Tools_test.nf'
include { DOWNLOAD_SNPEFF_DB } from './modules/references/Test/Tools_test.nf'
include { DOWNLOAD_ARRIBA } from './modules/references/Test/Tools_test.nf'
include { DOWNLOAD_VEP_CACHE } from './modules/references/Test/VEP_test.nf'
include { DOWNLOAD_CLINVAR } from './modules/references/Test/VEP_test.nf'


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
        referenceGenomePath = "${params.test_data_dir}/reference/genome.fa"
    } else {
        referenceGenomePath = genome_path.toString()
    }
    println " Reference genome path set to: ${referenceGenomePath}"
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
		genomeIndexPath = "${params.test_data_dir}/reference/genome.fa.fai"
	} else {
        genomeIndexPath = genome_index_path.toString()
    }
    println " Genome index path set to: ${genomeIndexPath}"
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
        genomeDictPath = "${params.test_data_dir}/reference/genome.dict"
    } else {
        genomeDictPath = genome_dict_path.toString()
    }
    println " Genome dictionary path set to: ${genomeDictPath}"
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
        gtfPath = "${params.test_data_dir}/reference/annotations.gtf"
    } else {
        gtfPath = gtf_path.toString()
    }
    println " GTF annotation path set to: ${gtfPath}"
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
        starIndexPath = "${params.test_data_dir}/reference/STAR_index"
    } else {
        starIndexPath = star_index_path.toString()
    }
    println " STAR genome index path set to: ${starIndexPath}"
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
        denylistPath = "${params.test_data_dir}/reference/denylist.bed"
    } else {
        denylistPath = denylist_path.toString()
    }
    println " Denylist BED file path set to: ${denylistPath}"
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
        snpVcfPath = "${params.test_data_dir}/reference/variants_snp.vcf.gz"
    } else {
        snpVcfPath = snp_vcf_path.toString()
    }
    println " SNP VCF path set to: ${snpVcfPath}"
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
        snpIndexPath = "${params.test_data_dir}/reference/variants_snp.vcf.gz.tbi"
    } else {
        snpIndexPath = snp_index_path.toString()
    }
    println "SNP Index path set to: ${snpIndexPath}"
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
        indelsVcfPath = "${params.test_data_dir}/reference/variants_indels.vcf.gz"
    } else {
        indelsVcfPath = indels_vcf_path.toString()
    }
    println " Indels VCF path set to: ${indelsVcfPath}"
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
        indelsIndexPath = "${params.test_data_dir}/reference/variants_indels.vcf.gz.tbi"
    } else {
        indelsIndexPath = indels_index_path.toString()
    }
    println " Indels Index path set to: ${indelsIndexPath}"
}


	// ========================== Filter and Merge VCFs ========================== //
    
	
	// Define channels for merged VCF and index
    def merged_vcf_ch
    def merged_vcf_index_ch

    if (file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz").exists() &&
        file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi").exists()) {

        //Creating channels if it already exists
        merged_vcf_ch = Channel.of(file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz"))
        merged_vcf_index_ch = Channel.of(file("${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"))

        println " Merged VCF and index already exist in the publish directory. Skipping merge."

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
        mergedVcfPath = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz"
    } else {
        mergedVcfPath = merged_vcf_path.toString()
    }
    println " Merged VCF path set to: ${mergedVcfPath}"
}

	// Capture merged VCF index path
	merged_vcf_index_ch.view { merged_vcf_index_path ->
    if (merged_vcf_index_path.toString().contains('/work/')) {
        mergedVcfIndexPath = "${params.test_data_dir}/reference/merged.filtered.recode.vcf.gz.tbi"
    } else {
        mergedVcfIndexPath = merged_vcf_index_path.toString()
    }
    println "Merged VCF Index path set to: ${mergedVcfIndexPath}"
}

	
	// ========================== Java Version Check ========================== //

    def java_check_ch

    if (file("${params.test_data_dir}/reference/java_check.log").exists()) {
        println " Java check log already exists in the publish directory. Skipping Java check."

        // If the log exists, create a channel from the existing log file
        java_check_ch = Channel.of(file("${params.test_data_dir}/reference/java_check.log"))

    } else {
        // Run the CHECK_JAVA process if the log file does not exist
        java_check_ch = CHECK_JAVA().java_output
    }

    // ========================== Capture Java Check Output ========================== //

    java_check_ch.view { java_log_path ->
        javaLogPath = file("${params.test_data_dir}/reference/java_check.log").toString()
        println "Java check log path set to: ${javaLogPath}"
    }
	
	// ========================== SnpEff Tool Handling ========================== //

    def snpeff_jar_ch, snpeff_config_ch

	if (file(params.snpeff_jar_path).exists() && file(params.snpeff_config_path).exists()) {
		println " SnpEff tool found in the server directory. Skipping download."

		snpeff_jar_ch = Channel.of(file(params.snpeff_jar_path))
		snpeff_config_ch = Channel.of(file(params.snpeff_config_path))

		snpEffJarPath = params.snpeff_jar_path
		snpEffConfigPath = params.snpeff_config_path

	} else if (file("${params.test_data_dir}/Tools/snpEff/snpEff.jar").exists() && 
			file("${params.test_data_dir}/Tools/snpEff/snpEff.config").exists()) {
    
		println " SnpEff tool found in the publish directory. Skipping download."

		snpeff_jar_ch = Channel.of(file("${params.test_data_dir}/Tools/snpEff/snpEff.jar"))
		snpeff_config_ch = Channel.of(file("${params.test_data_dir}/Tools/snpEff/snpEff.config"))

		snpEffJarPath = "${params.test_data_dir}/Tools/snpEff/snpEff.config"
		snpEffConfigPath = "${params.test_data_dir}/Tools/snpEff/snpEff.config"

	} else {
		println " SnpEff tool not found. Downloading..."
		def result = DOWNLOAD_SNPEFF_TOOL()
		snpeff_jar_ch = result.snpeff_jar
		snpeff_config_ch = result.snpeff_config
}


    // ========================== Capture SnpEff Paths ========================== //

    // Capture SnpEff JAR path
	snpeff_jar_ch.view { snpeff_jar_path ->  
    if (snpeff_jar_path.toString().contains('/work/')) {
        snpEffJarPath = "${params.test_data_dir}/Tools/snpEff/snpEff.jar"
    } else {
        snpEffJarPath = snpeff_jar_path.toString()
    }
    println " SnpEff JAR path set to: ${snpEffJarPath}"
}

	// Capture SnpEff Config path
	snpeff_config_ch.view { snpeff_config_path ->  
    if (snpeff_config_path.toString().contains('/work/')) {
        snpEffConfigPath = "${params.test_data_dir}/Tools/snpEff/snpEff.config"
    } else {
        snpEffConfigPath = snpeff_config_path.toString()
    }
    println "SnpEff Config path set to: ${snpEffConfigPath}"
}

	
	// ========================== SnpEff Database Handling ========================== //

	def snpeff_db_ch  

	if (params.snpeff_db_dir && file("${params.snpeff_db_dir_path}/${params.genomedb}").exists()) {
    println " SnpEff database for ${params.genomedb} found in the server directory. Skipping download."

    snpeff_db_ch = Channel.of(file("${params.snpeff_db_dir_path}/${params.genomedb}"))
    snpEffDbPath = "${params.snpeff_db_dir_path}/${params.genomedb}"

	} else if (file("${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}").exists()) {
    println " SnpEff database found in the publish directory. Skipping download."

    snpeff_db_ch = Channel.of(file("${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"))
    snpEffDbPath = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"

	} else {
    println " SnpEff database not found. Downloading..."
    def result = DOWNLOAD_SNPEFF_DB(params.genomedb, snpeff_jar_ch)

    snpeff_db_ch = result
    snpEffDbPath = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"
	}

// ========================== Capture SnpEff Database Path ========================== //

	// Capture SnpEff Database path
	snpeff_db_ch.view { snpeff_db_path ->  
    if (snpeff_db_path.toString().contains('/work/')) {
        snpEffDbPath = "${params.test_data_dir}/Tools/snpEff/snpEff/data/${params.genomedb}"
    } else {
        snpEffDbPath = snpeff_db_path.toString()
    }
    println " SnpEff Database path set to: ${snpEffDbPath}"
}



// ========================== Arriba Tool Handling ========================== //


	def arriba_dir_ch

	if (params.arriba_tool_dir_path && file("${params.arriba_tool_dir_path}/arriba_v2.4.0").exists()) {
    println " Arriba tool found in the server directory."
    arriba_dir_ch = Channel.of(file("${params.arriba_tool_dir_path}/arriba_v2.4.0"))

	} else if (file("${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0").exists()) {
    println " Arriba tool found in the publish directory."
    arriba_dir_ch = Channel.of(file("${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0"))

	} else {
    println " Arriba tool not found. Downloading..."
    def result = DOWNLOAD_ARRIBA()
    arriba_dir_ch = result.arriba_dir
	}

// ========================== Capture Arriba Tool and Database Paths ========================== //
// ========================== Handle Arriba Tool and Capture Paths ========================== //
	arriba_dir_ch.view { arriba_dir_path ->  
    // Set Arriba path after checking if it exists in the server or publish directory
    if (arriba_dir_path.toString().contains('/work/')) {
        arribaPath = "${params.test_data_dir}/Tools/ARRIBA/arriba_v2.4.0"
    } else {
        arribaPath = arriba_dir_path.toString()
    }

    println " Arriba tool path set to: ${arribaPath}"

    // Define known fusions and blacklist paths directly without immediate file check
    //knownFusionsPath = "${arribaPath}/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
    //blacklistPath = "${arribaPath}/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"

    //println " Known fusions path set to: ${knownFusionsPath}"
    //println " Blacklist path set to: ${blacklistPath}"
}




// ========================== VEP Cache Handling ========================== //
    def vep_cache_ch

    if (params.vep_cache_dir_path && file("${params.vep_cache_dir_path}").exists()) {
        println " VEP cache found in the server directory."
        vep_cache_ch = Channel.of(file("${params.vep_cache_dir_path}"))

    } else if (file("${params.test_data_dir}/Tools/VEP").exists()) {
        println " VEP cache found in the publish directory."
        vep_cache_ch = Channel.of(file("${params.test_data_dir}/Tools/VEP"))

    } else {
        println " VEP cache not found. Downloading..."
        def result = DOWNLOAD_VEP_CACHE()
        vep_cache_ch = result.vep_cache
    }

    // ========================== Capture VEP Cache Path ========================== //
	vep_cache_ch.view { vep_cache_path ->  
    if (vep_cache_path.toString().contains('/work/')) {
        vepCachePath = "${params.test_data_dir}/Tools/VEP"
    } else {
        vepCachePath = vep_cache_path.toString()
    }
    println " VEP Cache path set to: ${vepCachePath}"
}


// ========================== ClinVar VCF Handling ========================== //
    def clinvar_vcf_ch, clinvar_tbi_ch

    if (file("${params.clinvar_path}").exists() && file("${params.clinvartbi_path}").exists()) {
        println " ClinVar VCF and index found in the server directory."
        clinvar_vcf_ch = Channel.of(file("${params.clinvar_path}"))
        clinvar_tbi_ch = Channel.of(file("${params.clinvartbi_path}"))

    } else if (file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz").exists() && 
               file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi").exists()) {
        println " ClinVar VCF and index found in the publish directory."
        clinvar_vcf_ch = Channel.of(file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz"))
        clinvar_tbi_ch = Channel.of(file("${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi"))

    } else {
        println "️ ClinVar VCF and index not found. Downloading..."
        def result = DOWNLOAD_CLINVAR()
        clinvar_vcf_ch = result.out[0]
        clinvar_tbi_ch = result.out[1]
    }

    // ========================== Capture ClinVar Paths ========================== //
    clinvar_vcf_ch.view { clinvar_vcf_path ->  
    if (clinvar_vcf_path.toString().contains('/work/')) {
        clinvarVcfPath = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz"
    } else {
        clinvarVcfPath = clinvar_vcf_path.toString()
    }
    println " ClinVar VCF path set to: ${clinvarVcfPath}"
}

	clinvar_tbi_ch.view { clinvar_tbi_path ->  
    if (clinvar_tbi_path.toString().contains('/work/')) {
        clinvarTbiPath = "${params.test_data_dir}/Tools/VEP/clinvar.vcf.gz.tbi"
    } else {
        clinvarTbiPath = clinvar_tbi_path.toString()
    }
    println " ClinVar VCF Index path set to: ${clinvarTbiPath}"
}





    // ========================== Writing to Config After Completion ========================== //
    workflow.onComplete {
        try {
            def baseDir = System.getProperty('user.dir')
            def outputDir = new File("${baseDir}/reference_test_paths.config").getParentFile()

            if (!outputDir.exists()) {
                println " Output directory does not exist. Creating: ${outputDir}"
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
			params.vep_cache_dir = '${vepCachePath ?: 'NOT_FOUND'}'
			params.clinvar = '${clinvarVcfPath ?: 'NOT_FOUND'}'
            params.clinvartbi = '${clinvarTbiPath ?: 'NOT_FOUND'}'
            """

            println " Reference paths successfully written to ${configFile}"

        } catch (Exception e) {
            println " Error writing reference paths: ${e.message}"
        }
    }
}
