nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESSING } from './Subworkflows/preprocessing.nf'
include { VARIANT_CALLING } from './Subworkflows/variant_calling.nf'
include { ANNOTATION } from './Subworkflows/annotation.nf'
include { GENE_FUSION } from './Subworkflows/gene_fusion.nf'
include { MULTIQC_REPORT } from './Subworkflows/multiqc.nf'



workflow {

    if (params.only_variant_calling && params.only_fusion_detection) {
        error "You cannot enable both only_variant_calling and only_fusion_detection at the same time. Set one to false."
    }

    //==================================== Preprocessing ====================================//
    
    //==============================Step 1: Quality Control ===================================//
	
    PREPROCESSING(params.samplesheet)

    // Outputs from Preprocessing Subworkflow
	
    trimmed_reads_ch  = PREPROCESSING.out.trimmed_reads
    fastp_reports_ch  = PREPROCESSING.out.fastp_reports
    qc_results_ch     = PREPROCESSING.out.qc_reports
    multiqc_quality   = PREPROCESSING.out.multiqc
	

	
    if (params.only_qc) {
        log.info("Running only QC pipeline...")
        return
    }

    //==================================== Run Specific Pipelines ====================================//


        //===============================Step 2: Variant Calling=====================================//
		
    if (params.only_variant_calling) {
        log.info("Running only Variant Calling pipeline... ")

        VARIANT_CALLING(
            trimmed_reads_ch,
            params.star_genome_index,
            params.reference_genome,
            params.reference_genome_index,
            params.reference_genome_dict,
            params.merged_vcf,
            params.merged_vcf_index,
            params.denylist_bed
        )
		
		//=======================Outputs from Variant calling Subworkflow===========================//

        final_vcf_output = VARIANT_CALLING.out.final_vcf
        star_logs = VARIANT_CALLING.out.star_logs
        samtools_flagstat = VARIANT_CALLING.out.samtools_flagstat
        gatk_metrics = VARIANT_CALLING.out.gatk_metrics
        bcftools_stats = VARIANT_CALLING.out.bcftools_stats
		filtered_vcf_stats = VARIANT_CALLING.out.filtered_vcf_stats
		
        log.info " Variant Calling Pipeline Completed."
		
		//==============================Step 3: Annotation======================================//
		
		log.info("🔬 Running Variant Annotation pipeline...")

    ANNOTATION(
		VARIANT_CALLING.out.final_vcf,     
		params.merge_vcf,
		params.snpeff_jar,
		params.snpeff_config,
		params.snpeff_db,
		params.genomedb,
		params.vep_cache_dir,
		params.clinvar,
		params.clinvartbi
	)



    log.info " Variant Annotation Pipeline Completed."
	
	//===========================================Step 4: MultiQC===================================//
	
	// Run MultiQC ONLY for Variant Calling
        
		log.info(" Running MultiQC for Variant Calling...")
        
		MULTIQC_REPORT(
            fastqc_results = qc_results_ch ,
            fastp_reports = fastp_reports_ch,
            star_logs = star_logs,
            samtools_flagstat = samtools_flagstat,
            gatk_metrics = gatk_metrics,
            bcftools_stats = bcftools_stats,
			filtered_vcf_stats = filtered_vcf_stats
        )
        log.info "MultiQC Report Generated for Variant Calling."
		
    }
	
	//======================================== Step 5: Gene Fusion Detection===========================//
    
	else if (params.only_fusion_detection) {
        log.info(" Running only Gene Fusion Detection pipeline...")

        GENE_FUSION(
            trimmed_reads_ch,
            params.star_genome_index,
            params.reference_genome,
            params.gtf_annotation,
            params.arriba_blacklist,
            params.arriba_known_fusions,
            params.scripts_dir
        )

        log.info " Gene Fusion Pipeline Completed."
    }
	
	//=====================================  Step 6: Run Both (Full Workflow)=============================//
    else {
        log.info(" Running both Variant Calling & Gene Fusion pipelines...")

        // Run Variant Calling
        VARIANT_CALLING(
            trimmed_reads_ch,
            params.star_genome_index,
            params.reference_genome,
            params.reference_genome_index,
            params.reference_genome_dict,
            params.merged_vcf,
            params.merged_vcf_index,
            params.denylist_bed
        )

        // Run Annotation
        ANNOTATION(
            VARIANT_CALLING.out.final_vcf,     
            params.merge_vcf,
            params.snpeff_jar,
            params.snpeff_config,
            params.snpeff_db,
            params.genomedb,
            params.vep_cache_dir,
            params.clinvar,
            params.clinvartbi
        )
		
		//  Run MultiQC ONLY for Variant Calling
        log.info(" Running MultiQC for Variant Calling...")
        
        MULTIQC_REPORT(
            fastqc_results = qc_results_ch ,
            fastp_reports = fastp_reports_ch,
            star_logs = star_logs,
            samtools_flagstat = samtools_flagstat,
            gatk_metrics = gatk_metrics,
            bcftools_stats = bcftools_stats,
			filtered_vcf_stats = filtered_vcf_stats
        )

        log.info " MultiQC Report Generated for Variant Calling."

        // Run Gene Fusion Detection
        GENE_FUSION(
            trimmed_reads_ch,
            params.star_genome_index,
            params.reference_genome,
            params.gtf_annotation,
            params.arriba_blacklist,
            params.arriba_known_fusions,
            params.scripts_dir
        )

        log.info " Complete Pipeline Execution Completed."
		
	
    }

}

	
	
    
