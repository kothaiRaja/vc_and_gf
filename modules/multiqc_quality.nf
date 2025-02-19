process MultiQC_quality {
  publishDir "${params.outdir}/multiqc_quality", mode: "copy"
  container "https://depot.galaxyproject.org/singularity/multiqc%3A1.24.1--pyhdfd78af_0"
  input:
    path report_files
  output:
    path "multiqc_report.html", emit: report
  """
  multiqc ${report_files.join(' ')} -o .
  """
}