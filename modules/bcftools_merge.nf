//====================================================Merging all the vcf files============================================//

process BCFTOOLS_MERGE {
	container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
	publishDir "${params.outdir}/BCFTOOLS_MERGE", mode: "copy"
	
    input:
    path vcfs

    output:
    path "merged_output.vcf.gz"
    path "merged_output.vcf.gz.tbi"
	path "variants_per_sample.tsv" 


    script:
    """
    vcfs_absolute=\$(for vcf in ${vcfs.join(' ')}; do realpath "\$vcf"; done)

    echo "Using absolute paths:" > debug_absolute_paths.log
    echo \$vcfs_absolute >> debug_absolute_paths.log

    bcftools merge \\
        \$vcfs_absolute \\
        -O z -o merged_output.vcf.gz
    tabix -p vcf merged_output.vcf.gz
	
	# Query merged VCF to generate a tabular summary
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t]\\n' merged_output.vcf.gz > variants_per_sample.tsv
    """
}