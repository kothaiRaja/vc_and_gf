nextflow.enable.dsl = 2



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
	path filtered_vcf_stats

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
	cp ${filtered_vcf_stats} multiqc_input/ || true

    ls -lh multiqc_input/

    multiqc multiqc_input -o .
    """
}