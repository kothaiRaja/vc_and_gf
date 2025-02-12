nextflow.enable.dsl = 2

process CONCAT_FASTQ {
    tag { sample_id }
    container "https://depot.galaxyproject.org/singularity/ubuntu%3A24.04"
    publishDir "${params.outdir}/concatenated", mode: "copy"

    input:
    tuple val(sample_id), path(reads), val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_R1.merged.fastq.gz"), path("${sample_id}_R2.merged.fastq.gz"), val(strandedness), emit: merged_reads

    script:
    """
    echo "Processing reads for sample: ${sample_id}"
    ls -l ${reads.join(' ')}

    R1_FILES=(\$(ls ${reads.join(' ')} | grep -E '_R1|read1' || true))
    R2_FILES=(\$(ls ${reads.join(' ')} | grep -E '_R2|read2' || true))

    if [[ \${#R1_FILES[@]} -gt 1 || \${#R2_FILES[@]} -gt 1 ]]; then
        echo "Multiple files detected. Concatenating..."
        cat \${R1_FILES[@]} | gzip > ${sample_id}_R1.merged.fastq.gz
        cat \${R2_FILES[@]} | gzip > ${sample_id}_R2.merged.fastq.gz
        echo "Merging complete for sample: ${sample_id}"
    else
        echo "Only one file per read direction detected. Skipping concatenation..."
        cp \${R1_FILES[0]:-${reads[0]}} ${sample_id}_R1.merged.fastq.gz
        cp \${R2_FILES[0]:-${reads[1]}} ${sample_id}_R2.merged.fastq.gz
        echo "Files copied as-is for sample: ${sample_id}"
    fi
    """
}

process FASTQC_RAW {
    tag { sample_id }
    publishDir "${params.outdir}/fastqc", mode: "copy"
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
    container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.4--h125f33a_5"
    publishDir "${params.outdir}/fastp", mode: "copy"

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
      --length_required 10 \
      --cut_front --cut_tail \
      --cut_window_size 10 \
      --cut_mean_quality 10 \
      --html ${sample_id}_fastp.html \
      --json ${sample_id}_fastp.json
    """
}

process STAR_ALIGNMENT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/STAR/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(trimmed_r1), path(trimmed_r2), val(strandedness)
    path star_index_dir

    output:
    tuple val(sample_id), path("${sample_id}_Aligned.sortedByCoord.out.bam"), path("${sample_id}_Log.final.out"), 
          path("${sample_id}_Log.out"), 
          path("${sample_id}_Log.progress.out"),
          val(strandedness)

    script:
    """
    STAR --genomeDir $star_index_dir \
     --readFilesIn ${trimmed_r1} ${trimmed_r2} \
     --runThreadN $task.cpus \
     --readFilesCommand zcat \
     --outFilterType BySJout \
     --alignSJoverhangMin 5 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMatchNmin 16 \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterScoreMinOverLread 0.3 \
     --outSAMattributes NH HI AS nM MD NM \
     --outFileNamePrefix ${sample_id}_pass1_

    STAR --genomeDir $star_index_dir \
     --readFilesIn ${trimmed_r1} ${trimmed_r2} \
     --runThreadN $task.cpus \
     --readFilesCommand zcat \
     --sjdbFileChrStartEnd ${sample_id}_pass1_SJ.out.tab \
     --outFilterType BySJout \
     --alignSJoverhangMin 5 \
     --alignSJDBoverhangMin 1 \
     --outFilterMismatchNmax 999 \
     --outFilterMatchNmin 16 \
     --outFilterMatchNminOverLread 0.3 \
     --outFilterScoreMinOverLread 0.3 \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 16000000000 \
     --outSAMattrRGline ID:$sample_id LB:library PL:illumina PU:machine SM:$sample_id \
     --outSAMunmapped Within \
     --outSAMattributes NH HI AS nM MD NM \
     --outFileNamePrefix ${sample_id}_

    """
}

process SAMTOOLS_SORT_INDEX {
	tag { sample_id }
    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/sorted_bam", mode: "copy"

    input:
    tuple val(sample_id), path(bam), val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), path("${sample_id}_sorted.bam.bai"), val(strandedness)

    script:
    """
    # Sort the BAM file
    samtools sort -o ${sample_id}_sorted.bam ${bam}

    # Index the sorted BAM file
    samtools index ${sample_id}_sorted.bam
    """
}

process SAMTOOLS_FILTER_ORPHANS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--hd87286a_0"
    publishDir "${params.outdir}/filtered_bam", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai),val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.bam"), path("${sample_id}_filtered.bam.bai"),val(strandedness)

    script:
    """
    # Filter out orphan reads (retain properly paired reads)
    samtools view -h -f 2 ${sorted_bam} | samtools view -b - > ${sample_id}_filtered.bam

    # Index the filtered BAM file
    samtools index ${sample_id}_filtered.bam
    """
}

process SAMTOOLS_FLAGSTAT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"
    publishDir "${params.outdir}/flagstat", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(bai), val(strandedness)

    output:
    tuple val(sample_id), path("${sample_id}_flagstat.txt"), path("${sample_id}_stats_report.txt"), val(strandedness)

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

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/mark_duplicates", mode: "copy"

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bam_index),val(strandedness)

    output:
    tuple val(sample_id), 
          path("${sample_id}_marked_duplicates.bam"), 
          path("${sample_id}_marked_duplicates.bai"),
		  val(strandedness),
          path("${sample_id}_dup_metrics.txt")

    script:
    """
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


process SPLIT_NCIGAR_READS {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/split_ncigar", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai),val(strandedness)
	path(genome_fasta)
	path (index)
	path (genome_dict)

    output:
    tuple val(sample_id), 
          path("${sample_id}_split.bam"), 
          path("${sample_id}_split.bai"),
		  val(strandedness)

    script:
    """
    gatk SplitNCigarReads \
        -R ${genome_fasta} \
        -I ${bam} \
        -O ${sample_id}_split.bam \
        --create-output-bam-index true  
		
	"""
}



process SAMTOOLS_CALMD {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/samtools%3A1.18--h50ea8bc_1"

    publishDir "${params.outdir}/calmd", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai),val(strandedness) 
          path(genome_fasta)
		  path (index)
		  
		  
    output:
    tuple val(sample_id), 
          path("${sample_id}_calmd.bam"), 
          path("${sample_id}_calmd.bam.bai"),
		  val(strandedness)
         

    script:
    """
    # Add NM and MD tags using samtools calmd
    samtools calmd -b ${bam} ${genome_fasta} > ${sample_id}_calmd.bam

    # Index the updated BAM file
    samtools index ${sample_id}_calmd.bam

    """
}


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


process BED_TO_INTERVAL_LIST {
    tag { bed_file.name }

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/intervals", mode: "copy"

    input:
    path bed_file
    path genome_fasta
    path genome_dict

    output:
    path("${bed_file.baseName}.interval_list")

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
    path interval_list
    path genome_dict

    output:
    path "*.interval_list"

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


process GATK_HAPLOTYPE_CALLER {
    tag { sample_id }

     container "https://depot.galaxyproject.org/singularity/gatk4%3A4.4.0.0--py36hdfd78af_0"
    publishDir "${params.outdir}/haplotype_caller", mode: "copy"

    input:
    tuple val(sample_id), path(bam), path(bai) ,val(strandedness)
	path (genome)
	path (genome_index)
	path (genome_dict)
	path (interval_list)

    output:
    tuple val(sample_id), path("output_${sample_id}.vcf.gz"), path("output_${sample_id}.vcf.gz.tbi"), val(strandedness)

    script:
	def intervals_args = interval_list.collect { "--intervals ${it}" }.join(' ')
	
    """
    # Validate Inputs
    if [ ! -s ${bam} ]; then
        echo "Error: BAM file not found or empty." >&2
        exit 1
    fi
    if [ ! -s ${genome} ]; then
        echo "Error: Reference genome not found or empty." >&2
        exit 1
    fi

    # Run HaplotypeCaller
    gatk HaplotypeCaller \
        --native-pair-hmm-threads ${task.cpus} \
        --reference ${genome} \
        --output output_${sample_id}.vcf.gz \
        -I ${bam} \
        --standard-min-confidence-threshold-for-calling 30.0 \
        --dont-use-soft-clipped-bases true \
        --min-base-quality-score 10 \
        --output-mode EMIT_VARIANTS_ONLY \
        --read-filter OverclippedReadFilter \
        ${intervals_args} \
        --verbosity INFO
    """

  
}


process BCFTOOLS_STATS {
    tag "bcftools_stats"

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h8b25389_0" 
    publishDir "${params.outdir}/bcftools_stats", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index) , val(strandedness)  

    output:
    tuple val(sample_id), val(strandedness), path("haplotypecaller_stats_${sample_id}.txt")   

    script:
    """
    # Generate stats
    bcftools stats output_${sample_id}.vcf.gz > haplotypecaller_stats_${sample_id}.txt

    """
}

//=====================================================Individual vcf files=======================================================//
process GATK_VARIANT_FILTER {
    tag "${sample_id}_variant_filter"

    container "https://depot.galaxyproject.org/singularity/gatk4%3A4.2.6.0--hdfd78af_0"
    publishDir "${params.outdir}/variant_filter/${sample_id}", mode: "copy"

    input:
    tuple val(sample_id), path(vcf_file), path(vcf_index), val(strandedness)
    path genome
    path genome_index
    path genome_dict

    output:
    tuple val(sample_id), path("${sample_id}_filtered.vcf.gz"), path("${sample_id}_filtered.vcf.gz.tbi"), val(strandedness)


    script:
    """
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

    # Redirect logs
    if [ ! -s ${sample_id}_filtered.vcf.gz ]; then
        echo "Error: Filtered VCF is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}




process BCFTOOLS_QUERY {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/bcftools%3A1.19--h8b25389_1"
    publishDir "${params.outdir}/variant_summary", mode: "copy"

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_index), val(strandedness)

    output:
    tuple val(sample_id), val(strandedness), path("filtered_variants_summary_${sample_id}.txt")

    script:
    """
    # Generate a summary of filtered variants
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' ${filtered_vcf} > filtered_variants_summary_${sample_id}.txt

   # Handle cases where the output file is empty
    if [ ! -s filtered_variants_summary_${sample_id}.txt ]; then
        echo "No variants detected for ${sample_id}" > filtered_variants_summary_${sample_id}.txt
        exit 0
    fi
    """
}

process ANNOTATE_INDIVIDUAL_VARIANTS {
    tag "${sample_id}_annotate"

    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index), val(strandedness)
	path (snpEffJar)
	path (snpEffConfig)
	path (snpEffDbDir)
	val (genomedb)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf"),val(strandedness), path("${sample_id}.annotated.summary.html")

    script:
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        ${filtered_vcf} > ${sample_id}.annotated.vcf
		



    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats ${sample_id}.annotated.summary.html \
        ${filtered_vcf} > /dev/null
    """
}


process ANNOTATE_INDIVIDUAL_VARIANTS_VEP {
    tag "${sample_id}_vep_annotate"

    container "https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_index), val(strandedness) // Input VCF and its index
    path vep_cache                                                  // VEP cache directory
    path clinvar_vcf                                                // ClinVar VCF file
    path clinvar_index                                              // ClinVar Tabix index file

    output:
    tuple val(sample_id), path("${sample_id}.vep.annotated.vcf"),val(strandedness), path("${sample_id}.vep.summary.html")

    script:
    """
    # Annotate using Ensembl VEP with ClinVar
    vep --input_file ${filtered_vcf} \
        --output_file ${sample_id}.vep.annotated.vcf \
        --stats_file ${sample_id}.vep.summary.html \
        --cache \
        --dir_cache ${vep_cache} \
        --assembly GRCh38 \
        --format vcf \
        --vcf \
        --symbol \
        --protein \
        --force_overwrite \
        --custom ${clinvar_vcf},ClinVar,vcf,exact,0,CLNSIG,CLNDN

    # Validate that the annotated VCF is not empty
    if [ ! -s ${sample_id}.vep.annotated.vcf ]; then
        echo "Error: VEP annotation output is empty for ${sample_id}" >&2
        exit 1
    fi
    """
}

process EXTRACT_individual_VCF {

	tag { sample_id }
    container params.container
	publishDir "${params.outdir}/annotations", mode: 'copy'
	
	
    input:
    tuple val(sample_id), path(vcf_file),val(strandedness)

    output:
    path '*'

    script:
    """
    python - <<EOF
import pandas as pd
import vcfpy

# File paths
vcf_file_path = '${vcf_file}'
sample_id = '${sample_id}'
strandedness = '${strandedness}'

# Initialize list to store parsed VCF data
vcf_rows = []

# Parse the VCF file
reader = vcfpy.Reader.from_path(vcf_file_path)

# Extract headers for samples
sample_names = reader.header.samples.names

# Extract data from VCF
for record in reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = ",".join(str(a) for a in record.ALT)
    filt = ",".join(record.FILTER)
    dp = record.INFO.get('DP', 'NA')  # Extract Depth of Coverage
    ann_field = record.INFO.get('ANN', ['NA'])

    if ann_field != 'NA':
        ann_details = [a.split('|') for a in ann_field]
        for ann in ann_details:
            gene = ann[3] if len(ann) > 3 else "NA"
            impact = ann[2] if len(ann) > 2 else "NA"
            variant_type = ann[1] if len(ann) > 1 else "NA"

            # Iterate over samples to include sample-specific data
            for sample in sample_names:
                sample_data = record.call_for_sample[sample].data if record.call_for_sample.get(sample) else None
                genotype = sample_data.get('GT', './.') if sample_data else './.'

                vcf_rows.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "FILTER": filt,
                    "DP": dp,
                    "Impact": impact,
                    "Gene": gene,
                    "Variant": variant_type,
                    "Sample_ID": sample,
                    "Genotype": genotype,
					"Strandedness": strandedness 
                })

# Convert the VCF data to a DataFrame
vcf_df = pd.DataFrame(vcf_rows)

# Remove exact duplicate rows
vcf_df = vcf_df.drop_duplicates()

# Save VCF data to a file for reference
vcf_output_path = f"{sample_id}_parsed_annotations.csv"
vcf_df.to_csv(vcf_output_path, index=False)

print(f"Extracted VCF data saved to: {vcf_output_path}")

EOF
"""
    
}




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


process ANNOTATE_VARIANTS {
    tag "Annotate variants"
    
    container "https://depot.galaxyproject.org/singularity/snpeff%3A5.2--hdfd78af_1"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf	// Filtered VCF file
	path index
	path tsv
    path snpEffJar  // Path to the SnpEff JAR file
    path snpEffConfig  // Path to the SnpEff configuration file
    path snpEffDbDir	// Path to the SnpEff database directory
	val genomedb

    output:
    path "annotated.vcf"
	path "annotated.summary.html"

    script:
	
    
    """
    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        ${vcf} > annotated.vcf

    java -Xmx16G -jar ${snpEffJar} \
        -c ${snpEffConfig} \
        -v ${params.genomedb} \
        -dataDir ${snpEffDbDir} \
        -stats annotated.summary.html \
        ${vcf} > /dev/null
    """
}


process EXTRACT_VCF {
    container params.container
	publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf_file
	path html

    output:
    path "extracted_vcf_data.csv"

    script:
    """
    python - <<EOF
import pandas as pd
import vcfpy

# File paths
vcf_file_path = "${vcf_file}"

# Initialize list to store parsed VCF data
vcf_rows = []

# Parse the VCF file
reader = vcfpy.Reader.from_path(vcf_file_path)

# Extract headers for samples
sample_names = reader.header.samples.names

# Extract data from VCF
for record in reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = ",".join(str(a) for a in record.ALT)
    filt = ",".join(record.FILTER)
    dp = record.INFO.get('DP', 'NA')  # Extract Depth of Coverage
    ann_field = record.INFO.get('ANN', ['NA'])

    if ann_field != 'NA':
        ann_details = [a.split('|') for a in ann_field]
        for ann in ann_details:
            gene = ann[3] if len(ann) > 3 else "NA"
            impact = ann[2] if len(ann) > 2 else "NA"
            variant_type = ann[1] if len(ann) > 1 else "NA"

            # Iterate over samples to include sample-specific data
            for sample in sample_names:
                sample_data = record.call_for_sample[sample].data
                genotype = sample_data.get('GT', './.') if sample_data else './.'

                vcf_rows.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "FILTER": filt,
                    "DP": dp,
                    "Impact": impact,
                    "Gene": gene,
                    "Variant": variant_type,
                    "Sample_ID": sample,
                    "Genotype": genotype
                })

# Convert the VCF data to a DataFrame
vcf_df = pd.DataFrame(vcf_rows)

# Save VCF data to a file for reference
vcf_output_path = "extracted_vcf_data.csv"
vcf_df.to_csv(vcf_output_path, index=False)

print(f"Extracted VCF data saved to: {vcf_output_path}")
EOF
    """
}


process ANNOTATEVARIANTS_VEP {
    container 'https://depot.galaxyproject.org/singularity/ensembl-vep%3A110.1--pl5321h2a3209d_0'
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path input_vcf          // Input VCF file
    path input_vcf_tbi      // Tabix index file for the VCF
	path tsv
    path vep_cache          // VEP cache directory
    path clinvar_vcf        // ClinVar VCF file
    path clinvar_vcf_tbi    // Tabix index file for ClinVar
	

    output:
    path "annotated_variants.vcf"          // Annotated VCF output
    path "annotated_variants.html"         // HTML summary report

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
        --input_file ${input_vcf} \
        --output_file annotated_variants.vcf \
        --stats_file annotated_variants.html \
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

process MULTIQC_REPORT {
    tag "MultiQC Report"

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
    [ -e "${fastqc_results}" ] && cp ${fastqc_results} multiqc_input/
    [ -e "${fastp_reports}" ] && cp ${fastp_reports} multiqc_input/
    [ -e "${star_logs}" ] && cp ${star_logs} multiqc_input/
    [ -e "${samtools_flagstat}" ] && cp ${samtools_flagstat} multiqc_input/
    [ -e "${gatk_metrics}" ] && cp ${gatk_metrics} multiqc_input/
    [ -e "${bcftools_stats}" ] && cp ${bcftools_stats} multiqc_input/

    # Run MultiQC
    multiqc multiqc_input -o .
    """
}

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
    samples_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample_id, [file(row.fastq_1), file(row.fastq_2)], row.strandedness) }

    concatenated_reads_ch = params.concatenate ? 
        CONCAT_FASTQ(samples_ch) : 
        samples_ch.map { sample_id, reads, strandedness -> tuple(sample_id, reads[0], reads[1], strandedness) }

    qc_results_ch = FASTQC_RAW(concatenated_reads_ch)
    trimmed_reads_ch = TRIM_READS(concatenated_reads_ch)
	
	// Extract FastQC file paths (flatten both .zip and .html files)
	qc_files_ch = qc_results_ch.map { it[1] + it[2] }.flatten()

	// Extract Fastp file paths (exclude sample name)
	fastp_files_ch = trimmed_reads_ch.fastp_reports.map { it[1..-1] }.flatten()

	// Combine FastQC and Fastp channels
	combined_channel = qc_files_ch.concat(fastp_files_ch).collect()

	// Pass combined file paths to MultiQC
	multiqc_quality = MultiQC_quality(combined_channel)

	trimmed_reads_ch.fastp_reports.view { "Fastp reports passed to MultiQC: ${it}" }


    if (params.only_qc) {
        qc_results_ch.view { "QC results for sample: ${it}" }
    } else{
        star_aligned_ch = STAR_ALIGNMENT(trimmed_reads_ch[0], params.star_genome_index)
		//Step 3: Sort and index BAM files
		sorted_bams = SAMTOOLS_SORT_INDEX(star_aligned_ch.map { tuple(it[0], it[1], it[5]) })

		//Step 4: Filter the Orphan reads
		filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams)
	
		//Step 5: Generate alignment statistics
		alignment_stats = SAMTOOLS_FLAGSTAT(filtered_bams)
	
		//Step 6: Mark duplicates
		marked_bams = GATK_MARK_DUPLICATES(filtered_bams)
		
		marked_bams.view { "Mark Duplicates Output: $it" }

	
		// Step 7: Split N CIGAR reads
		split_bams = SPLIT_NCIGAR_READS(marked_bams.map { 
		tuple(it[0], it[1], it[2], it[3]) }, params.reference_genome, params.reference_genome_index, params.reference_genome_dict)

		// Step 8: SAMTOOLS CALMD process
		calmd_ch = SAMTOOLS_CALMD(split_bams.map { tuple(it[0], it[1], it[2], it[3]) }, params.reference_genome, params.reference_genome_index) 

		// Step 9: Recalibrate and Apply BQSR
		recalibrated_bams = GATK_RECALIBRATION(calmd_ch.map { 
		tuple(it[0], it[1], it[2], it[3]) }, params.reference_genome, params.reference_genome_index, params.reference_genome_dict, params.merged_vcf, params.merged_vcf_index)

		
		//Step 10: Convert BED to interval list
		interval_list_ch = BED_TO_INTERVAL_LIST(params.denylist_bed, params.reference_genome, params.reference_genome_dict)
	
		//Step 11: Scatter the Interval List
		scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, params.reference_genome_dict)
	
		//Step 12: GATK HaplotypeCaller
		gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams.map{tuple(it[0], it[1], it[2], it[3]) }, params.reference_genome, params.reference_genome_index, params.reference_genome_dict, scattered_intervals_ch)
				
	
		gvcf_output.view { "Raw GVCF output: $it" }
	
		//Step 13: provide stats
		bcftools_stats_ch = BCFTOOLS_STATS(gvcf_output)
		
		// Step 15: Filter individual VCF files
		filtered_individual_vcfs = GATK_VARIANT_FILTER(gvcf_output, params.reference_genome, params.reference_genome_index, params.reference_genome_dict)
		
		// Step 16: Provide Stats
		filtered_vcf_stats = BCFTOOLS_QUERY(filtered_individual_vcfs)
		
		// Step 17: Create a mapped collection of filtered VCF paths for merging
		filtered_vcf = filtered_individual_vcfs
					.map { it[1] } // Extract paths to filtered VCF files
					.collect()     // Collect them into a list

		filtered_vcf.view { "Filtered VCF output: $it" }
		
		// Step 19: Conditional processing for merging or annotating individual files
		if (params.merge_vcf) {
		// Merge filtered VCFs
		merged_filtered_vcfs = BCFTOOLS_MERGE(filtered_vcf)
		
		// Step 20: Annotate the merged VCF file snpeff
		annotated_merged_vcf = ANNOTATE_VARIANTS(merged_filtered_vcfs, params.snpeff_jar, params.snpeff_config, params.snpeff_db, params.genomedb)
	
		// Step 21: Annotate the merged vcfs ensembl_vep
		annotated_merged_vcf_ensemblvep = ANNOTATEVARIANTS_VEP(merged_filtered_vcfs, params.vep_cache_dir, params.clinvar, params.clinvartbi)

		println "Merging and annotating VCF files completed."
		
		// Step 22: Create a table from the annotated merged VCF
		table_creation = EXTRACT_VCF(annotated_merged_vcf)

		println "Table creation from merged VCF completed."
		
		} else {
    
	// Step 23: Annotate individual VCF files
    annotated_individual_vcfs = ANNOTATE_INDIVIDUAL_VARIANTS(filtered_individual_vcfs, params.snpeff_jar, params.snpeff_config, params.snpeff_db, params.genomedb)
	
	// Step 24: Annotate individual VCF files ensemblvep
	annotated_individual_vcf_ensemblvep = ANNOTATE_INDIVIDUAL_VARIANTS_VEP (filtered_individual_vcfs, params.vep_cache_dir, params.clinvar, params.clinvartbi)
	
	//Step25: individual_vcf to table
	csv_individual_vcf = EXTRACT_individual_VCF(annotated_individual_vcfs.map {tuple(it[0], it[1], it[2])})

    println "Individual VCF annotation completed."
    println "Table creation step skipped because merging is disabled."
}

	multiqc_results = MULTIQC_REPORT(
    fastqc_results = qc_results_ch.map { [it[1], it[2]] }.flatten(),   
    fastp_reports = trimmed_reads_ch.fastp_reports.map { [it[1], it[2]] }.flatten(),
    star_logs = star_aligned_ch.map { [it[2], it[3], it[4]] }.flatten(),
    samtools_flagstat = alignment_stats.map { it[1] },  
    gatk_metrics = marked_bams.map { it[4] },  
    bcftools_stats = bcftools_stats_ch.map { it[2] } 
	)
     
		println "MultiQC report generated: ${multiqc_results}"
	
	//==================================================Gene Fusion===============================================//



//Step 27: STAR alignment for fusion detection
	star_align_fusion_ch = STAR_ALIGN_FUSION(trimmed_reads_ch[0], params.star_genome_index, params.gtf_annotation )

//Step 28: Fusion detection using ARRIBA
	ARRIBA_ch = ARRIBA(star_align_fusion_ch, params.reference_genome, params.gtf_annotation, params.known_blacklist_chr22, params.known_fusion_chr22 )
	
//Step 29: Visualization step
    fusion_visuals = ARRIBA_VISUALIZATION(ARRIBA_ch, params.scripts_dir, params.reference_genome, params.gtf_annotation)	
}

}