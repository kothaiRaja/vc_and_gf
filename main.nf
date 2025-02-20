nextflow.enable.dsl = 2

// Import subworkflows
include { PREPROCESSING } from './Subworkflows/preprocessing.nf'
include { VARIANT_CALLING } from './Subworkflows/variant_calling.nf'
include { ANNOTATION } from './Subworkflows/annotation.nf'


process MULTIQC_REPORT {
    container "https://depot.galaxyproject.org/singularity/multiqc%3A1.14--pyhdfd78af_0"
    publishDir "${params.outdir}/multiqc", mode: "copy"

    input:
    path fastqc_results
    path fastp_reports
    path star_logs
    path samtools_flagstat
    path gatk_metrics
    path bcftools_stats

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    mkdir -p multiqc_input

    # Copy available input files if they exist
    cp ${fastqc_results} multiqc_input/ || true
    cp ${fastp_reports} multiqc_input/ || true
    cp ${star_logs} multiqc_input/ || true
    cp ${samtools_flagstat} multiqc_input/ || true
    cp ${gatk_metrics} multiqc_input/ || true
    cp ${bcftools_stats} multiqc_input/ || true

    ls -lh multiqc_input/

    multiqc multiqc_input -o .
    """
}





//===========================================STAR Fusion=============================================//
process STAR_ALIGN_FUSION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
	publishDir "${params.outdir}/Star_fusion", mode: 'copy'
    
	input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
    path star_index_dir            
    path gtf                                          

    output:
    output:
    tuple val(sample_id), path('*Log.final.out'), val(strandedness), emit: log_final
    tuple val(sample_id), path('*Log.out'), val(strandedness), emit: log_out
    tuple val(sample_id), path('*Aligned.sortedByCoord.out.bam'), val(strandedness), emit: bam_sorted
    tuple val(sample_id), path('*Chimeric.out.sam'), val(strandedness), emit: chimeric_sam
    tuple val(sample_id), path('*Log.progress.out'), val(strandedness), emit: log_progress
    tuple val(sample_id), path('*SJ.out.tab'), val(strandedness), emit: splice_junctions

    
    script:
    """
    STAR --genomeDir $star_index_dir \
     --readFilesIn ${trimmed_r1} ${trimmed_r2} \
     --runThreadN 8 \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMatchNmin 16 \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterScoreMinOverLread 0.3 \
     --chimSegmentMin 10 \
     --chimJunctionOverhangMin 10 \
	 --chimScoreJunctionNonGTAG -4 \
	 --chimScoreMin 1 \
     --chimOutType WithinBAM SeparateSAMold \
     --chimScoreDropMax 50 \
     --chimScoreSeparation 10 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM MD NM \
	 --limitBAMsortRAM 32000000000 \
     --outFileNamePrefix ${sample_id}_

    """
}


process ARRIBA {
    tag { sample_id }
    
    container "https://depot.galaxyproject.org/singularity/arriba%3A2.4.0--hdbdd923_3"
    publishDir "${params.outdir}/ARRIBA", mode: 'copy'

    input:
	tuple val(sample_id),path(log_final), val(strandedness)
	tuple val(sample_id),path(log_out), val(strandedness)
	tuple val(sample_id),path(bam), val(strandedness)
	tuple val(sample_id),path(chimeric_sam), val(strandedness)
	tuple val(sample_id), path(log_progress) , val(strandedness)
    tuple val(sample_id), path(splice_junctions), val(strandedness) 
    path fasta                                     
    path gtf                                       
    path blacklist                                 
    path known_fusions  
                           

    output:
	tuple val(sample_id), path("${sample_id}.fusions.tsv"), val(strandedness), emit: fusions
	tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), val(strandedness), emit: fusions_discarded


    script:
    """
    arriba \
         -x $bam \
         -c ${sample_id}_Chimeric.out.sam \
         -a $fasta \
         -g $gtf \
         -b $blacklist \
         -k $known_fusions \
         -o ${sample_id}.fusions.tsv \
         -O ${sample_id}.fusions.discarded.tsv
    """
}


process ARRIBA_VISUALIZATION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/r-base%3A4.4.1"  
    publishDir "${params.outdir}/ARRIBA_VISUALIZATION", mode: 'copy'

    input:
    tuple val(sample_id), path(fusions_tsv), val(strandedness)
	tuple val (sample_id), path(discarded_tsv), val(strandedness)
    path r_script  
    path fasta  
    path gtf
  

    output:
    path "*.fusion_plot.pdf"

    script:
    """
    # Extract prefix (filename without extension)
    PREFIX=\$(basename ${fusions_tsv} .fusions.tsv)

    # Debug: Log resolved paths
    echo "Resolved fusions file path: \$(realpath ${fusions_tsv})"
    echo "Resolved annotation file path: \$(realpath ${gtf})"
    

    # Verify the resolved paths exist
    if [ ! -f "\$(realpath ${fusions_tsv})" ]; then
        echo "Error: Fusions file not found at: \$(realpath ${fusions_tsv})"
        exit 1
    fi
    if [ ! -f "\$(realpath ${gtf})" ]; then
        echo "Error: Annotation file not found at: \$(realpath ${gtf})"
        exit 1
    fi
	
	# Check if the TSV file has more than just headers
    DATA_LINES=\$(awk 'NR > 1 {print; exit}' ${fusions_tsv})
    if [ -z "\$DATA_LINES" ]; then
        echo "No fusions detected for ${sample_id}. Skipping visualization." >&2
        touch \${PREFIX}.fusion_plot.pdf
        exit 0
    fi

    # Run the R script
    Rscript draw_fusions.R \
        --fusions "\$(realpath ${fusions_tsv})" \
        --annotation "\$(realpath ${gtf})" \
        --output "\${PREFIX}.fusion_plot.pdf"
    
    """
}






workflow {

    if (params.only_variant_calling && params.only_fusion_detection) {
        error "You cannot enable both only_variant_calling and only_fusion_detection at the same time. Set one to false."
    }

    //==================================== Preprocessing ====================================//
    
    // ✅ Correctly pass `samplesheet` into the PREPROCESSING workflow
    PREPROCESSING(params.samplesheet)

    // ✅ Assign outputs correctly
    trimmed_reads_ch  = PREPROCESSING.out.trimmed_reads
    fastp_reports_ch  = PREPROCESSING.out.fastp_reports
    qc_results_ch     = PREPROCESSING.out.qc_reports
    multiqc_quality   = PREPROCESSING.out.multiqc
	

	
    if (params.only_qc) {
        log.info("Running only QC pipeline...")
        return
    }

    //==================================== Run Specific Pipelines ====================================//


        // **Step 2: Variant Calling**
    if (params.only_variant_calling) {
        log.info("Running only Variant Calling pipeline... 🧬")

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

        final_vcf_output = VARIANT_CALLING.out.final_vcf
        stats_report = VARIANT_CALLING.out.stats_report

        log.info "✅ Variant Calling Pipeline Completed."
		
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



    log.info "✅ Variant Annotation Pipeline Completed."
    }
}

	
	
    
