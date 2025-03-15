nextflow.enable.dsl = 2

process FASTQC_RAW {
    tag { sample_id }
	
	cpus params.get('fastqc_raw_cpus', 2)
    memory params.get('fastqc_raw_memory', '4 GB')
    time params.get('fastqc_raw_time', '1h')
	
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_fastqc.*"
    container "https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0"

    input:
    tuple val(sample_id), path(r1), path(r2), val(strandedness)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html"), val(strandedness), emit: qc_results

    script:
    """
    fastqc ${r1} ${r2} --outdir .
    """
}

process TRIM_READS {
    tag { sample_id }
	
	cpus params.get('trim_reads_cpus', 4)
	memory params.get('trim_reads_memory', '8 GB')
	time params.get('trim_reads_time', '2h')

	
    container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--h125f33a_5"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_fastp.*"

    input:
    tuple val(sample_id), path(r1), path(r2), val(strandedness)

    output:
    tuple val(sample_id), path("trimmed_${sample_id}_R1.fastq.gz"), path("trimmed_${sample_id}_R2.fastq.gz"), val(strandedness), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp.html"), path("${sample_id}_fastp.json"), emit: fastp_reports

    script:
    """
    fastp -i ${r1} -I ${r2} \
      -o trimmed_${sample_id}_R1.fastq.gz \
      -O trimmed_${sample_id}_R2.fastq.gz \
      --detect_adapter_for_pe \
      --adapter_sequence auto \
      --adapter_sequence_r2 auto \
      --length_required 50 \
      --cut_front --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --html ${sample_id}_fastp.html \
      --json ${sample_id}_fastp.json
    """
}


process STAR_ALIGNMENT {
    tag { sample_id }

    cpus params.get('star_cpus', 12)
    memory params.get('star_memory', '24GB')
	time params.get('star_time', '6h')

    container params.get('star_container', "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0")
    publishDir "${params.outdir}/star_output", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
    path star_index_dir
    path gtf_file

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_id), val(strandedness), path("${sample_id}_Log.final.out"), emit: log_final
    tuple val(sample_id), val(strandedness), path("${sample_id}_Log.out"), emit: log_out
    tuple val(sample_id), val(strandedness), path("${sample_id}_Log.progress.out"), emit: log_progress
    tuple val(sample_id), val(strandedness), path("${sample_id}_Chimeric.out.sam"), optional: true, emit: chimeric_sam
    tuple val(sample_id), val(strandedness), path("${sample_id}_SJ.out.tab"), optional: true, emit: junctions


    script:
    def extra_args = params.get('star_extra_args', '') 
    def out_sam_type = "BAM SortedByCoordinate"
	def out_sam_attr = "--outSAMattrRGline ID:${sample_id} LB:library PL:${params.get('seq_platform', 'ILLUMINA')} PU:machine SM:${sample_id} CN:${params.get('seq_center', 'Unknown')}"
	def RAM_LIMIT = task.memory.toMega() * 1000000
	// Determine if we should use --outSAMstrandField intronMotif
    def strand_option = ""
    if (strandedness == 'unstranded') {
        strand_option = "--outSAMstrandField intronMotif"
    } 


    """
    echo "Running STAR Alignment for Sample: ${sample_id}"

    THREADS=${task.cpus}
    RAM_LIMIT=${RAM_LIMIT}



    # **First Pass: Discover Splice Junctions**
    STAR --genomeDir ${star_index_dir} \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN \$THREADS \
         --readFilesCommand zcat \
         --outFilterType BySJout \
         --alignSJoverhangMin ${params.get('star_alignSJoverhangMin', 8)} \
         --alignSJDBoverhangMin ${params.get('star_alignSJDBoverhangMin', 1)} \
         --outFilterMismatchNmax ${params.get('star_outFilterMismatchNmax', 999)} \
         --outFilterMatchNmin ${params.get('star_outFilterMatchNmin', 16)} \
         --outFilterMatchNminOverLread ${params.get('star_outFilterMatchNminOverLread', 0.3)} \
         --outFilterScoreMinOverLread ${params.get('star_outFilterScoreMinOverLread', 0.3)} \
         --outSAMattributes NH HI AS nM MD NM \
         --outFileNamePrefix ${sample_id}_pass1_ \
         --limitBAMsortRAM \$RAM_LIMIT \
		 $strand_option \
         

    # **Second Pass: Final Alignment**
    STAR --genomeDir ${star_index_dir} \
         --readFilesIn ${trimmed_r1} ${trimmed_r2} \
         --runThreadN \$THREADS \
         --readFilesCommand zcat \
         --sjdbFileChrStartEnd ${sample_id}_pass1_SJ.out.tab \
         --sjdbGTFfile ${gtf_file} \
         --twopassMode None \
         --outFilterType BySJout \
         --alignSJoverhangMin ${params.get('star_alignSJoverhangMin', 8)} \
         --alignSJDBoverhangMin ${params.get('star_alignSJDBoverhangMin', 1)} \
         --outFilterMismatchNmax ${params.get('star_outFilterMismatchNmax', 999)} \
         --outFilterMatchNmin ${params.get('star_outFilterMatchNmin', 16)} \
         --outFilterMatchNminOverLread ${params.get('star_outFilterMatchNminOverLread', 0.3)} \
         --outFilterScoreMinOverLread ${params.get('star_outFilterScoreMinOverLread', 0.3)} \
         --chimSegmentMin ${params.get('star_chimSegmentMin', 10)} \
         --chimJunctionOverhangMin ${params.get('star_chimJunctionOverhangMin', 10)} \
         --chimScoreJunctionNonGTAG ${params.get('star_chimScoreJunctionNonGTAG', -4)} \
         --chimScoreMin ${params.get('star_chimScoreMin', 1)} \
         --chimOutType ${params.get('star_chimOutType', 'WithinBAM SeparateSAMold')} \
         --chimScoreDropMax ${params.get('star_chimScoreDropMax', 50)} \
         --chimScoreSeparation ${params.get('star_chimScoreSeparation', 10)} \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM \$RAM_LIMIT \
         --outSAMunmapped Within \
		 --quantMode TranscriptomeSAM GeneCounts \
		 --outFileNamePrefix ${sample_id}_ \
         ${out_sam_attr} \
		 $strand_option \
		$extra_args
    """
}




process SAMTOOLS_SORT_INDEX {
	tag { sample_id }
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/sorted_bam", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai")

    script:
    """
    # Sort the BAM file
    samtools sort -o ${sample_id}_sorted.bam ${bam}

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam
    """
}

process SAMTOOLS_STATS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_stats_report.*"

    input:
    tuple val(sample_id), val(strandedness), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_stats_report.txt")
	
    script:
    """
    # Run samtools stats on the sorted BAM file
    samtools stats ${sorted_bam} > ${sample_id}_stats_report.txt
    """
}


process SAMTOOLS_FILTER_ORPHANS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/filtered_bam", mode: "copy"

    input:
    tuple val(sample_id),val(strandedness), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id),val(strandedness), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai")

    script:
    """
    # Filter out orphan reads (retain properly paired reads)
    samtools view -h -F 3844 ${sorted_bam} | samtools view -b - > ${sample_id}_filtered.bam


	# Index the filtered BAM file
    samtools index ${sample_id}_filtered.bam
    """
}



process SAMTOOLS_FLAGSTAT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_flagstat.*"

    input:
    tuple val(sample_id), val(strandedness), path(sorted_bam), path(bai)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_flagstat.txt"), path("${sample_id}_stats_report.txt")

    script:
    """
    # Validate BAM and index file
    if [ ! -s ${sorted_bam} ]; then
        echo "Error: BAM file is empty or missing for ${sample_id}" >&2
        exit 1
    fi
    
    if [ ! -s ${bai} ]; then
        echo "Error: BAM index (.bai) file is missing for ${sample_id}" >&2
        exit 1
    fi

    # Validate BAM file to ensure it is not corrupted
    samtools quickcheck -v ${sorted_bam} || (echo "BAM file validation failed for ${sample_id}" && exit 1)

    # Run samtools flagstat to generate alignment metrics
    samtools flagstat ${sorted_bam} > ${sample_id}_flagstat.txt

    # Generate additional quality metrics with samtools stats
    samtools stats ${sorted_bam} > ${sample_id}_stats_report.txt
    """
}

process GATK_MARK_DUPLICATES {
    tag { sample_id }
	
	cpus params.get('gatk_mark_duplicates_cpus', 8)
    memory params.get('gatk_mark_duplicates_memory', '16 GB')
    time params.get('gatk_mark_duplicates_time', '2h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/dedup_bam", mode: "copy"

    input:
    tuple val(sample_id),val(strandedness), path(sorted_bam), path(sorted_bam_index)

    output:
     output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_marked_duplicates.bam"), path("${sample_id}_marked_duplicates.bai")
	tuple val(sample_id), val(strandedness), path("${sample_id}_dup_metrics.txt")

    script:
    """
    THREADS=${task.cpus}

    gatk MarkDuplicates \
        -I ${sorted_bam} \
        -O ${sample_id}_marked_duplicates.bam \
        -M ${sample_id}_dup_metrics.txt \
        --CREATE_INDEX true \
		--REMOVE_DUPLICATES ${params.remove_duplicates ? 'true' : 'false'} \
        --VALIDATION_STRINGENCY ${params.validation_stringency ?: 'LENIENT'}
		
	# Check if output BAM is generated
    if [ ! -s ${sample_id}_marked_duplicates.bam ]; then
        echo "Error: Marked duplicates BAM file not generated for ${sample_id}" >&2
        exit 1
    fi

    # Check if metrics file is generated
    if [ ! -s ${sample_id}_dup_metrics.txt ]; then
        echo "Error: Duplicate metrics file not generated for ${sample_id}" >&2
        exit 1
    fi
    """
}

process BED_TO_INTERVAL_LIST {
    
    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/intervals", mode: "copy"

    input:
    tuple val(meta_id), path(bed_file)
	path(genome_fasta)
	path(genome_dict)

    output:
    tuple val(meta_id), path("${bed_file.baseName}.interval_list")

    script:
    """
    gatk BedToIntervalList \
        -I ${bed_file} \
        -O ${bed_file.baseName}.interval_list \
        -SD ${genome_dict}
    """
}

process SCATTER_INTERVAL_LIST {
    tag "Scatter interval list"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/scattered_intervals", mode: 'copy'

    input:
    tuple val(meta), path(interval_list)
	path(genome_dict)

    output:
    tuple val(meta), path("*.interval_list")

    script:
    """
    mkdir -p scattered_intervals
    gatk IntervalListTools \
        --INPUT ${interval_list} \
        --OUTPUT scattered_intervals \
        --SCATTER_COUNT ${params.scatter_count} \
        --UNIQUE true

    # Move and rename the output files to the working directory with unique names
    for f in scattered_intervals/*/*; do
        dir_name=\$(basename \$(dirname "\$f"))  # Get subdirectory name
        file_name=\$(basename "\$f")            # Get original file name
        mv "\$f" "scattered_intervals/\${dir_name}_\${file_name}.interval_list"
    done

    # Move the uniquely renamed files to the working directory
    mv scattered_intervals/*.interval_list .
    rm -r scattered_intervals
	
	# Log the generated intervals
    echo "Scattered intervals:"
    cat *.interval_list
    """
}

process SPLIT_NCIGAR_READS {
    tag { "${sample_id}_${interval.baseName}" }  

    cpus params.get('split_ncigar_reads_cpus', 8)
    memory params.get('split_ncigar_reads_memory', '16 GB')
    time params.get('split_ncigar_reads_time', '3h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"

    publishDir "${params.outdir}/split_ncigar", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
    path genome_fasta
    path index
    path genome_dict

    output:
    tuple val(sample_id),val(strandedness), 
          path("${sample_id}_split_${interval.baseName}.bam"), 
          path("${sample_id}_split_${interval.baseName}.bai")
          

    script:
    def avail_mem = task.memory ? task.memory.giga : 3  
    def interval_command = interval ? "--intervals ${interval}" : ""
    

    """
    echo "Running SplitNCigarReads for sample: ${sample_id} on interval: ${interval.baseName}"

    gatk --java-options "-Xmx${avail_mem}g" SplitNCigarReads \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        -O ${sample_id}_split_${interval.baseName}.bam \\
        --skip-mapping-quality-transform false \\
        --max-mismatches-in-overhang 1 \\
        --max-bases-in-overhang 50 \\
        --create-output-bam-index true \\
        --process-secondary-alignments true \\
        ${interval_command} 

    echo "SplitNCigarReads completed for sample: ${sample_id} on interval: ${interval.baseName}"

    
    """
}

process MERGE_BAMS {
    tag "MERGE_BAMS"

    cpus params.get('merge_bam_cpus', 4)
    memory params.get('merge_bam_memory', '8 GB')
    time params.get('merge_bam_time', '2h')

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.15.1--h1170115_0"

    input:
    tuple val(sample_id), val(strandedness), path(bam_list), path(bai_list)

    output:
    tuple val(sample_id),val(strandedness), 
          path("${sample_id}_merged.bam"), 
          path("${sample_id}_merged.bam.bai")
          

    script:
    """
    echo "Merging BAM files for sample: ${sample_id}"

    # Create a temporary file to store BAM paths
    ls ${bam_list} > bam_files.txt

    # Ensure the file contains BAM files
    if [ ! -s bam_files.txt ]; then
        echo "ERROR: No BAM files found for sample: ${sample_id}"
        exit 1
    fi

    # Merge BAM files
    samtools merge -@ ${task.cpus} -b bam_files.txt -o ${sample_id}_merged.bam

    # Index BAM file
    samtools index ${sample_id}_merged.bam

    echo "Merge complete for sample: ${sample_id}"
    """
}


process SAMTOOLS_CALMD {
    tag { "${sample_id}_${bam.baseName}" }  // Ensure unique tags for each split BAM

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"

    publishDir "${params.outdir}/calmd", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai) 
    path genome_fasta
    path index

    output:
    tuple val(sample_id),val(strandedness), 
          path("${bam.baseName}_calmd.bam"),
          path("${bam.baseName}_calmd.bam.bai")
          

    script:
    """
    echo "Processing BAM file: ${bam} with samtools calmd"

    # Add NM and MD tags using samtools calmd
    samtools calmd -b ${bam} ${genome_fasta} > ${bam.baseName}_calmd.bam

    # Verify BAM file exists before proceeding
    if [ ! -s ${bam.baseName}_calmd.bam ]; then
        echo "Error: samtools calmd did not produce a valid BAM file!" >&2
        exit 1
    fi

    # Ensure BAM is sorted before indexing
    samtools sort -o ${bam.baseName}_calmd.sorted.bam ${bam.baseName}_calmd.bam
    mv ${bam.baseName}_calmd.sorted.bam ${bam.baseName}_calmd.bam

    # Index the BAM file
    samtools index ${bam.baseName}_calmd.bam

    # Verify BAI file exists
    if [ ! -s ${bam.baseName}_calmd.bam.bai ]; then
        echo "Error: samtools index did not generate a BAI file!" >&2
        exit 1
    fi

    echo "Finished processing ${bam.baseName}_calmd.bam"
    """
}

process SAMTOOLS_IDXSTATS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_idxstats.*"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_idxstats.txt")
	
    script:
    """
    samtools idxstats ${bam} > ${sample_id}_idxstats.txt
    """
}


process GATK_BASERECALIBRATOR {
    tag { sample_id }
	
    cpus params.get('gatk_recalibration_cpus', 8)
    memory params.get('gatk_recalibration_memory', '32 GB')
    time params.get('gatk_recalibration_time', '5h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/recal_tables", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(interval)
    path(genome_fasta)
    path(index)
    path(dict)
    path(known_variants)
    path(known_variants_index)

    output:
    tuple val(sample_id), val(strandedness), path("${sample_id}_recal_data.table")

    script:

    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=${task.cpus}

    # Step 1: BaseRecalibrator
    gatk BaseRecalibrator \
        -R ${genome_fasta} \
        -I ${bam} \
        --known-sites ${known_variants} \
        -O ${sample_id}_recal_data.table \
		${interval_command} 
		
       

    # Check if recalibration table is generated
    if [ ! -s ${sample_id}_recal_data.table ]; then
        echo "Error: Recalibration table not generated for ${sample_id}" >&2
        exit 1
    fi
    """
}

process GATK_APPLYBQSR {
    tag { sample_id }

    cpus params.get('gatk_applybqsr_cpus', 8)
    memory params.get('gatk_applybqsr_memory', '32 GB')
    time params.get('gatk_applybqsr_time', '5h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/recalibrated_bams", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai), path(recal_table), path(interval)
    path(genome_fasta)
    path(index)
    path(dict)

    output:
    tuple val(sample_id),val(strandedness),	
          path("${sample_id}_${interval.baseName}_recalibrated.bam"),
          path("${sample_id}_${interval.baseName}_recalibrated.bai")
          

    script:
    
    def interval_command = interval ? "--intervals ${interval}" : ""

    """
    THREADS=${task.cpus}

    # Step 2: ApplyBQSR
    gatk ApplyBQSR \\
        -R ${genome_fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${sample_id}_${interval.baseName}_recalibrated.bam \\
        $interval_command

    # Check if recalibrated BAM is generated
    if [ ! -s ${sample_id}_${interval.baseName}_recalibrated.bam ]; then
        echo "Error: Recalibrated BAM not generated for ${sample_id}" >&2
        exit 1
    fi

    # Step 3: Index BAM using GATK
    gatk BuildBamIndex \\
        -I ${sample_id}_${interval.baseName}_recalibrated.bam
    """
}


process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }
    
    cpus params.get('gatk_haplotype_caller_cpus', 12)
    memory params.get('gatk_haplotype_caller_memory', '24 GB')
    time params.get('gatk_haplotype_caller_time', '8h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), val(strandedness), path(bam), path(bai)
    path(genome)
    path(genome_index)
    path(genome_dict)
    path(known_sites_vcf)
    path(known_sites_vcf_index)

    output:
    tuple val(sample_id), val(strandedness), path("output_${bam.baseName}.vcf.gz"), 
          path("output_${bam.baseName}.vcf.gz.tbi")

    script:
    """
    THREADS=${task.cpus}

    # Validate Inputs
    if [ ! -s ${bam} ]; then
        echo "Error: BAM file not found or empty." >&2
        exit 1
    fi
    if [ ! -s ${genome} ]; then
        echo "Error: Reference genome not found or empty." >&2
        exit 1
    fi

    # Extract a unique identifier from BAM filename
    BAM_BASENAME=\$(basename ${bam} .bam)

    # Run HaplotypeCaller with RNA-seq optimizations
    gatk HaplotypeCaller \
        --native-pair-hmm-threads \$THREADS \
        --reference ${genome} \
        --output output_\${BAM_BASENAME}.vcf.gz \
        -I ${bam} \
        --standard-min-confidence-threshold-for-calling 10.0 \
        --min-base-quality-score 10 \
        --output-mode EMIT_VARIANTS_ONLY \
        --dont-use-soft-clipped-bases true \
        --disable-read-filter NotDuplicateReadFilter \
        --dbsnp ${known_sites_vcf} \
        --verbosity INFO
    """
}




process GATK_MERGEVCFS {
    tag { sample_id }

    cpus params.get('gatk_merge_vcfs_cpus', 4)
    memory params.get('gatk_merge_vcfs_memory', '8GB')
    time params.get('gatk_merge_vcfs_time', '4h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/merged_vcf", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_list), path(tbi_list)

    output:
    tuple val(sample_id), path("merged_${sample_id}.vcf.gz"), path("merged_${sample_id}.vcf.gz.tbi")

    script:
    """
    gatk MergeVcfs \\
        ${vcf_list.collect { "-I ${it}" }.join(" \\\n")} \\
        -O merged_${sample_id}.vcf.gz

    gatk IndexFeatureFile -I merged_${sample_id}.vcf.gz
    """
}

process BCFTOOLS_STATS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.3--h577a1d6_9"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_bcftools_stats.*"

    input:
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}_bcftools_stats.txt")

    script:
    """
    bcftools stats ${vcf} > ${sample_id}_bcftools_stats.txt
    """
}

process GATK_VARIANT_FILTER {
    tag "${sample_id}_variant_filter"
    
    cpus params.get('gatk_variant_filter_cpus', 8)
    memory params.get('gatk_variant_filter_memory', '24 GB')
    time params.get('gatk_variant_filter_time', '3h')

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.6.1.0--py310hdfd78af_0"
    publishDir "${params.outdir}/variant_filter/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi")

    script:
    """
    THREADS=${task.cpus}

    # Run GATK VariantFiltration with configurable thresholds
    gatk VariantFiltration \
        -R ${genome} \
        -V ${vcf_file} \
        --cluster-window-size ${params.gatk_vf_window_size} \
        --cluster-size ${params.gatk_vf_cluster_size} \
        --filter-name "LowQual" --filter-expression "QUAL < ${params.gatk_vf_qual_filter}" \
        --filter-name "LowQD" --filter-expression "QD < ${params.gatk_vf_qd_filter}" \
        --filter-name "HighFS" --filter-expression "FS > ${params.gatk_vf_fs_filter}" \
        --filter-name "LowMQ" --filter-expression "MQ < ${params.gatk_vf_mq_filter}" \
        --filter-name "HighSOR" --filter-expression "SOR > ${params.gatk_vf_sor_filter}" \
        --filter-name "LowReadPosRankSum" --filter-expression "ReadPosRankSum < ${params.gatk_vf_read_pos_filter}" \
        --filter-name "LowBaseQRankSum" --filter-expression "BaseQRankSum < ${params.gatk_vf_baseq_filter}" \
        -O ${sample_id}_filtered.vcf.gz

    # Index the filtered VCF file for downstream compatibility
    gatk IndexFeatureFile -I ${sample_id}_filtered.vcf.gz

    # Validate output
    if [ ! -s ${sample_id}_filtered.vcf.gz ] || [ ! -s ${sample_id}_filtered.vcf.gz.tbi ]; then
        echo "Error: Filtered VCF or index is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}

process BCFTOOLS_QUERY {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.3--h577a1d6_9"
    publishDir "${params.outdir}/multiqc_input", mode: "copy", pattern: "*_query.*"

    input:
    tuple val(sample_id), path(filtered_vcf), path(tbi)

    output:
    tuple val(sample_id), path("${sample_id}_variant_summary.txt")

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n' ${filtered_vcf} > ${sample_id}_variant_summary.txt
    """
}

process SELECT_SNPs {
    tag "${sample_id}_select_snps"
	
	container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/selected_variants/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
	path index
	path dict

    output:
    tuple val(sample_id), path("${sample_id}_snps.vcf.gz"), path("${sample_id}_snps.vcf.gz.tbi")

    script:
    """
    gatk SelectVariants \
        -R ${genome} \
        -V ${vcf_file} \
        --select-type-to-include SNP \
        -O ${sample_id}_snps.vcf.gz

    gatk IndexFeatureFile -I ${sample_id}_snps.vcf.gz
    """
}

process SELECT_INDELs {
    tag "${sample_id}_select_indels"
	
	container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/selected_variants/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index)
    path genome
	path index
	path dict

    output:
    tuple val(sample_id), path("${sample_id}_indels.vcf.gz"), path("${sample_id}_indels.vcf.gz.tbi")

    script:
    """
    gatk SelectVariants \
        -R ${genome} \
        -V ${vcf_file} \
        --select-type-to-include INDEL \
        -O ${sample_id}_indels.vcf.gz

    gatk IndexFeatureFile -I ${sample_id}_indels.vcf.gz
    """
}


process ANNOTATE_VARIANTS {
    tag "Annotate variants"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path(snpEffJar)
    path(snpEffConfig)
    path(snpEffDbDir)
    val(genomedb)

    output:
    tuple val(sample_id), path("annotated_${sample_id}.vcf"), emit: annotated_vcf
    path "annotated_${sample_id}.summary.html", emit: summary_html

    script:
    """
    THREADS=${task.cpus}

    # Annotate variants and generate annotated VCF file
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated_${sample_id}.vcf

    # Generate annotation summary report
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}.summary.html \
        ${vcf} > /dev/null
    """
}

process ANNOTATE_VARIANTS_SELECTED {
    tag "${sample_id}_annotate"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotated_variants/${sample_id}", mode: "copy"

    cpus params.get('annotate_cpus', 8)
    memory params.get('annotate_memory', '16GB')
    time params.get('annotate_time', '2h')

    input:
    tuple val(sample_id), path(snp_vcf), path(snp_vcf_index), path(indel_vcf), path(indel_vcf_index)
    path(snpEffJar)
    path(snpEffConfig)
    path(snpEffDbDir)
    val(genomedb)

    output:
    tuple val(sample_id), path("${sample_id}_annotated_snps.vcf")
    tuple val(sample_id), path("${sample_id}_annotated_indels.vcf")
	path "annotated_${sample_id}_annotated_snps.summary.html"
	path "annotated_${sample_id}_annotated_indels.summary.html"
	

    script:
    """
	
	THREADS=${task.cpus}

    # Annotate variants and generate annotated VCF file
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${snp_vcf} > ${sample_id}_annotated_snps.vcf

   # Annotate variants and generate annotated VCF file
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${indel_vcf} > ${sample_id}_annotated_indels.vcf
		
	# Generate annotation summary report
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}_annotated_indels.summary.html \
        ${indel_vcf} > /dev/null
    
	# Generate annotation summary report
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated_${sample_id}_annotated_snps.summary.html \
        ${snp_vcf} > /dev/null
    
    """
}


process ANNOTATEVARIANTS_VEP {
    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_tbi)
    path vep_cache
	path clinvar_vcf        
    path clinvar_vcf_tbi
    

    output:
    tuple val(sample_id), path("annotated_${sample_id}.vcf"), emit: annotated_vcf  
    path "annotated_${sample_id}.html", emit: summary_html

    script:
    """
    # Debugging: Log resolved paths and existence
    echo "Input VCF: \$(realpath ${input_vcf})"
    echo "Input VCF Index: \$(realpath ${input_vcf_tbi})"
    echo "VEP Cache Directory: \$(realpath ${vep_cache})"
    echo "ClinVar VCF: \$(realpath ${clinvar_vcf})"
    echo "ClinVar Index: \$(realpath ${clinvar_vcf_tbi})"

    if [ ! -d "${vep_cache}" ]; then
        echo "Error: VEP cache directory does not exist at: \$(realpath ${vep_cache})" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf}" ]; then
        echo "Error: ClinVar VCF file not found at: \$(realpath ${clinvar_vcf})" >&2
        exit 1
    fi

    if [ ! -f "${clinvar_vcf_tbi}" ]; then
        echo "Error: ClinVar index file not found at: \$(realpath ${clinvar_vcf_tbi})" >&2
        exit 1
    fi

    # Run VEP annotation
    vep \
        --input_file "${input_vcf}" \
        --output_file "annotated_${sample_id}.vcf" \
        --stats_file "annotated_${sample_id}.html" \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN
    """

}

process EXTRACT_VARIANT_IMPACT {
    tag { sample_id }
    publishDir "${params.outdir}/annotations", mode: 'copy'
    container params.container

    cpus 2
    memory '2GB'
    time '30m'

    input:
    tuple val(sample_id), path(annotated_vcf)

    output:
    tuple val(sample_id), path("${sample_id}_variants.csv")

    script:
    """
    echo "Running Variant Extraction..."
   
    python3 <<EOF
import pandas as pd
import vcfpy

# File paths
vcf_file = "${annotated_vcf}"
output_csv = "${sample_id}_variants.csv"

# Initialize list to store parsed VCF data
variant_data = []

# Parse the VCF file
reader = vcfpy.Reader.from_path(vcf_file)

# Extract sample names
sample_names = reader.header.samples.names

# Extract data from VCF
for record in reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt_list = [str(a) for a in record.ALT] if record.ALT else ["NA"]
    filt = ",".join(record.FILTER) if record.FILTER else "PASS"
    dp = record.INFO.get('DP', 'NA')  

    # Determine variant type (SNV, Insertion, Deletion)
    for alt in alt_list:
        if len(ref) == len(alt) == 1:
            variant_type = "SNV"
        elif len(alt) > len(ref):
            variant_type = "Insertion"
        elif len(ref) > len(alt):
            variant_type = "Deletion"
        else:
            variant_type = "Complex"

        # Extract ANN field (annotation)
        ann_field = record.INFO.get('ANN', [])
        if ann_field:
            for annotation in ann_field:
                ann_details = annotation.split('|')
                gene = ann_details[3] if len(ann_details) > 3 else "NA"
                impact = ann_details[2] if len(ann_details) > 2 else "NA"
                effect = ann_details[1] if len(ann_details) > 1 else "NA"

                # Iterate over samples to include sample-specific data
                for sample in sample_names:
                    sample_data = record.call_for_sample.get(sample, None)
                    genotype = sample_data.data.get('GT', './.') if sample_data else './.'

                    variant_data.append({
                        "Sample_ID": sample,
                        "Chromosome": chrom,
                        "Position": pos,
                        "Reference": ref,
                        "Alternate": alt,
                        "Variant_Type": variant_type,
                        "Filter": filt,
                        "DP": dp,
                        "Impact": impact,
                        "Effect": effect,
                        "Gene": gene,
                        "Genotype": genotype
                    })

# Convert the VCF data to a DataFrame
vcf_df = pd.DataFrame(variant_data)

# Remove exact duplicate rows
vcf_df = vcf_df.drop_duplicates()

# Save to CSV
vcf_df.to_csv(output_csv, index=False)

print(f"Extracted variant data saved to: {output_csv}")
EOF
    """
}


process ARRIBA {
    tag { sample_id }

    cpus params.get('star_align_fusion_cpus', 12)
    memory params.get('star_align_fusion_memory', '24 GB')
    time params.get('star_align_fusion_time', '5h')

    container "https://depot.galaxyproject.org/singularity/arriba%3A2.4.0--hdbdd923_3"
    publishDir "${params.outdir}/ARRIBA", mode: 'copy'

    input:
    tuple val(sample_id), path(log_final), path(log_out), path(bam), path(chimeric_sam), path(log_progress), path(splice_junctions), val(strandedness)
    path fasta                                     
    path gtf                                       
    path blacklist                                 
    path known_fusions

    output:
    tuple val(sample_id), path("${sample_id}.fusions.tsv"), val(strandedness), emit: fusions
    tuple val(sample_id), path("${sample_id}.fusions.discarded.tsv"), val(strandedness), emit: fusions_discarded

    script:
    """
    echo "Running Arriba for Sample: ${sample_id}"

    arriba \
         -x $bam \
         -c $chimeric_sam \
         -a $fasta \
         -g $gtf \
         -b $blacklist \
         -k $known_fusions \
         -o ${sample_id}.fusions.tsv \
         -O ${sample_id}.fusions.discarded.tsv

    echo "Arriba finished for Sample: ${sample_id}"
    """
}


process ARRIBA_VISUALIZATION {
    tag { sample_id }
	
	cpus params.get('arriba_visualization_cpus', 4)
    memory params.get('arriba_visualization_memory', '16 GB')
    time params.get('arriba_visualization_time', '2h')


    container "https://depot.galaxyproject.org/singularity/r-base%3A4.4.1"  
    publishDir "${params.resultsdir}/ARRIBA_VISUALIZATION", mode: 'copy'

    input:
    tuple val(sample_id), path(fusions_tsv), val(strandedness)
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

// Step 1: Read samplesheet
    samples_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, file(row.fastq_1), file(row.fastq_2), row.strandedness) }


    // Step 3: Perform FastQC on raw reads
    qc_results_ch = FASTQC_RAW(samples_ch)

    // Step 4: Trim reads using Fastp
    trimmed_reads_ch = TRIM_READS(samples_ch)


	// **Step 1: STAR Alignment**
        star_aligned_ch = STAR_ALIGNMENT(trimmed_reads_ch[0], params.star_genome_index, params.gtf_annotation)
		
		


	

	
	

    // **Step 5: Sort BAM files**
    sorted_bams = SAMTOOLS_SORT_INDEX(star_aligned_ch[0])
    sorted_bams.view { "Sorted BAMs Output: $it" }
	

		//Step 3: samtools stats
		alignment_stats_ch = SAMTOOLS_STATS(sorted_bams)
		
	

        // **Step 3: Filter orphan reads**
        filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams)
		
		filtered_bams.view { "Filtered BAMs Output: $it" }
		
		

        // **Step 4: Generate alignment statistics**
        alignment_flagstats = SAMTOOLS_FLAGSTAT(filtered_bams)
		
		alignment_flagstats.view { "Alignment flagstats Output: $it" }
		

        // **Step 5: Mark duplicates using GATK**
        marked_bams = GATK_MARK_DUPLICATES(filtered_bams)
        
		
		
		//Step 6: Exons_BED channels
		ch_exon_bed = Channel
			.fromPath(params.exons_BED)
			.map { file_path -> tuple([id: file_path.baseName], file_path) }
			.view { meta, file -> println " ✅ Exons BED Tuple: ${meta}, File: ${file}" }

		
		
		
		// **Step 7: Convert BED to Interval List**
		interval_list_ch = BED_TO_INTERVAL_LIST(ch_exon_bed, params.reference_genome, params.reference_genome_dict)
		
		

		// **Step 10: Check if Scatter Intervals is Enabled**	
		scattered_intervals_ch = Channel.empty() 

		if (params.scatterintervals) {
    // Scatter the interval list
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.reference_genome_dict)
        .map { meta, file -> file }  // Remove tuple, keep only file paths
        .flatten()
    
    println "🔹 Using **scattered intervals** for parallel processing"
    
    // **VIEW statement to debug output**
    scattered_intervals_ch.view { file -> "Generated Scattered Interval: ${file}" }

	} else {
    // Use the full interval list directly
    scattered_intervals_ch = interval_list_ch.map { meta, file -> file }  // Remove tuple, keep only file paths
    println "🔹 Using **full exons interval list** without scattering"
	}

		// **Step 6: Pair Split BAMs with Scattered Intervals using combine()**
split_bams_ch = marked_bams[0]
    .map { tuple(it[0], it[1], it[2], it[3]) }
    .combine(scattered_intervals_ch)
    .map { sample_id, strandedness, bam, bai, interval -> 
        tuple(sample_id, strandedness, bam, bai, interval)
    }
    .set { ch_splitncigar_bam_bai_interval }

// Run SPLIT_NCIGAR_READS process using the paired inputs
split_bams = SPLIT_NCIGAR_READS(
    ch_splitncigar_bam_bai_interval,
    params.reference_genome, 
    params.reference_genome_index, 
    params.reference_genome_dict
)

split_bams.view { "SplitNcigar Output : $it" }

split_bams
    .groupTuple() // Groups by sample ID
    .map { sample_id, strandedness_list, bams, bais -> 
        tuple(sample_id, strandedness_list.unique()[0], bams.flatten(), bais.flatten()) 
    }
    .set { ch_merged_bams }


merged_bams = MERGE_BAMS(ch_merged_bams)



merged_bams.view { "Merged BAMs Output : $it" }



// **Step 7: SAMTOOLS CALMD**
        calmd_bams = SAMTOOLS_CALMD(merged_bams, params.reference_genome, params.reference_genome_index)
		
		calmd_bams.view { "Calmd BAM Output: $it" }
		
		// Now add idxstats here
idxstats_calmd_bams = SAMTOOLS_IDXSTATS(filtered_bams)
idxstats_calmd_bams.view { "Idxstats Output (filtered_bams): $it" }

	
	// Pair split BAMs with Scattered Intervals using combine()

	calmd_bams.map { tuple(it[0], it[1], it[2], it[3]) }
			.combine(scattered_intervals_ch)
			.map { sample_id, strandedness, bam, bai, interval -> 
			tuple(sample_id, strandedness, bam, bai, interval)
       
    }
	.set { ch_recalibrated_bam_bai_interval}
	
		//Run GATK Recalibration 
	recalibrated_bams_table = GATK_BASERECALIBRATOR(ch_recalibrated_bam_bai_interval,params.reference_genome, params.reference_genome_index, params.reference_genome_dict, params.merged_vcf, params.merged_vcf_index )
       
	recalibrated_bams_table.view{"Recalibrated Bams Output: $it" }
	
		
	
//Applying BSQR 
	
	// Join Calmd BAMs with Recalibrated Tables (Removing Redundant Strandedness)
ch_applybqsr = calmd_bams.map { sample_id, strandedness, calmd_bam, calmd_bai ->  
    tuple(sample_id, strandedness, calmd_bam, calmd_bai) 
}.join(recalibrated_bams_table.map { sample_id, _ , recal_table ->  // Ignore strandedness
    tuple(sample_id, recal_table)
}, by: 0)  // Join by sample_id

// Debug View to check joined data
ch_applybqsr.view { "Joined output: $it "}

// Combine with Scattered Intervals
ch_applybqsr.combine(scattered_intervals_ch)
    .map { sample_id, strandedness, bam, bai, recal_table, interval -> 
        tuple(sample_id, strandedness, bam, bai, recal_table, interval)
    }
    .set { ch_applybqsr_bam_bai_interval }

// Debug View to check final channel
ch_applybqsr_bam_bai_interval.view { "Joined combined channel: $it"}


//Run BSQR
	
	// Step 1: Apply BQSR on scattered BAMs
bams_base_recalibrated = GATK_APPLYBQSR(
    ch_applybqsr_bam_bai_interval, 
    params.reference_genome, 
    params.reference_genome_index, 
    params.reference_genome_dict 
)

bams_base_recalibrated.view { "recalibrated_bsqr_bams:$it " }


ch_GATK_variant_caller = GATK_HAPLOTYPE_CALLER(
    bams_base_recalibrated,
    params.reference_genome, 
    params.reference_genome_index, 
    params.reference_genome_dict, 
    params.merged_vcf, 
	params.merged_vcf_index
)

ch_GATK_variant_caller.view { "HaplotypeCaller_output: $it " }



ch_GATK_variant_caller
    .groupTuple() // Groups by sample ID
    .map { sample_id, strandedness_list, vcf, tbi -> 
        tuple(sample_id, strandedness_list.unique()[0], vcf.flatten(), tbi.flatten()) 
    }
    .set { ch_haplotypecaller_grouped }

// Debug Output
ch_haplotypecaller_grouped.view { "Grouped VCFs for Merging: $it" }


merged_vcf_ch = GATK_MERGEVCFS(
    ch_haplotypecaller_grouped.map { sample_id, strandedness, vcf_list, tbi_list -> 
        tuple(sample_id, vcf_list, tbi_list)
    }
)

bcftools_stats_ch = BCFTOOLS_STATS(merged_vcf_ch)


// Correct assignment of the new channel
ch_merged_vcfs = merged_vcf_ch.map { sample_id, vcf, tbi -> 
    tuple(sample_id, vcf, tbi)
}

// View the output
ch_merged_vcfs.view { it -> "Merged_vcf_output: $it" }


ch_variant_filter = GATK_VARIANT_FILTER(ch_merged_vcfs, params.reference_genome, params.reference_genome_index, params.reference_genome_dict)
        
ch_variant_filter.view {"Filtered_vcf: $it"}

bcftools_query_ch = BCFTOOLS_QUERY(ch_variant_filter)


if (params.select_variants) {
    selected_snps_ch = SELECT_SNPs(ch_variant_filter, params.reference_genome, params.reference_genome_index, params.reference_genome_dict)
	selected_snps_ch.view {"SNPs output: $it"}
    selected_indels_ch = SELECT_INDELs(ch_variant_filter, params.reference_genome, params.reference_genome_index, params.reference_genome_dict)
	selected_indels_ch.view {"INDELS output: $it"}

    // Combine SNPs & INDELs before annotation
        ch_selected_variants = selected_snps_ch.join(selected_indels_ch, by: 0)
		
		// View to debug the joined output
    ch_selected_variants.view { "Joined SNPs & INDELs: $it" }


    // Annotate both in the same process
    ch_annotated_variants = ANNOTATE_VARIANTS_SELECTED(ch_selected_variants,params.snpeff_jar, params.snpeff_config, params.snpeff_db, params.genomedb)
} else {
    ch_annotated_variants = ANNOTATE_VARIANTS(ch_variant_filter, params.snpeff_jar, params.snpeff_config, params.snpeff_db, params.genomedb)
	
	


    
    // Run the annotation process
    ANNOTATEVARIANTS_VEP(
        ch_variant_filter,
        params.vep_cache_dir, params.clinvar, params.clinvartbi  
        
    )
    


// **Run Extraction Process**
EXTRACT_VARIANT_IMPACT(ch_annotated_variants[0])
}



}