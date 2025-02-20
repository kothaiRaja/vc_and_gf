process GATK_RECALIBRATION {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/recalibrated_bams", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai),val(strandedness)
	path(genome_fasta)
	path(index)
	path(dict)
	path(known_variants)
	path(known_variants_index)

    output:
    tuple val(sample_id), 
          path("${sample_id}_recalibrated.bam"), 
          path("${sample_id}_recalibrated.bai"),
		  val(strandedness),
          path("${sample_id}_recal_data.table")
		  

    script:
    """
    # Step 1: BaseRecalibrator
    gatk BaseRecalibrator \
        -R ${genome_fasta} \
        -I ${bam} \
        --known-sites ${known_variants} \
        -O ${sample_id}_recal_data.table 
		
	# Check if recalibration table is generated
    if [ ! -s ${sample_id}_recal_data.table ]; then
        echo "Error: Recalibration table not generated for ${sample_id}" >&2
        exit 1
    fi
		

    # Step 2: ApplyBQSR
    gatk ApplyBQSR \
        -R ${genome_fasta} \
        -I ${bam} \
        --bqsr-recal-file ${sample_id}_recal_data.table \
        -O ${sample_id}_recalibrated.bam
		
	# Check if recalibrated BAM is generated
    if [ ! -s ${sample_id}_recalibrated.bam ]; then
        echo "Error: Recalibrated BAM not generated for ${sample_id}" >&2
        exit 1
    fi
	
	# Step 3: Index BAM using GATK
    gatk BuildBamIndex \
        -I ${sample_id}_recalibrated.bam

	# Step 4: Validation
	gatk ValidateSamFile \
		-I ${sample_id}_recalibrated.bam \
		-MODE SUMMARY > ${sample_id}_validation_summary.txt


    """
}
