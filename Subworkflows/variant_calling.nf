nextflow.enable.dsl = 2

// ✅ Include required modules
include { STAR_ALIGNMENT } from '../modules/star_alignment.nf'
include { SAMTOOLS_SORT_INDEX } from '../modules/samtools_sort_index.nf'
include { SAMTOOLS_FILTER_ORPHANS } from '../modules/samtools_filter_orphans.nf'
include { SAMTOOLS_FLAGSTAT } from '../modules/samtools_flagstat.nf'
include { GATK_MARK_DUPLICATES } from '../modules/gatk_mark_duplicates.nf'
include { SPLIT_NCIGAR_READS } from '../modules/split_ncigar_reads.nf'
include { SAMTOOLS_CALMD } from '../modules/samtools_calmd.nf'
include { GATK_RECALIBRATION } from '../modules/gatk_recalibration.nf'
include { BED_TO_INTERVAL_LIST } from '../modules/bed_to_interval_list.nf'
include { SCATTER_INTERVAL_LIST } from '../modules/scatter_interval_list.nf'
include { GATK_HAPLOTYPE_CALLER } from '../modules/gatk_haplotype_caller.nf'
include { GATK_VARIANT_FILTER } from '../modules/gatk_variant_filter.nf'
include { BCFTOOLS_STATS } from '../modules/bcftools_stats.nf'
include { BCFTOOLS_MERGE } from '../modules/bcftools_merge.nf'

workflow VARIANT_CALLING {
    take:
    trimmed_reads_ch
    star_index
    reference_genome
    genome_index
    genome_dict
    merged_vcf
    merged_vcf_index
    denylist_bed

    main:
    
    log.info "🧬 Starting Variant Calling Workflow..."

    // **Step 1: STAR Alignment**
    star_aligned_ch = STAR_ALIGNMENT(trimmed_reads_ch, star_index)

    // **Step 2: Sort BAM files**
    sorted_bams = SAMTOOLS_SORT_INDEX(star_aligned_ch.map { tuple(it[0], it[1], it[5]) })

    // **Step 3: Filter orphan reads**
    filtered_bams = SAMTOOLS_FILTER_ORPHANS(sorted_bams)

    // **Step 4: Generate alignment statistics**
    alignment_stats = SAMTOOLS_FLAGSTAT(filtered_bams)

    // **Step 5: Mark duplicates using GATK**
    marked_bams = GATK_MARK_DUPLICATES(filtered_bams)
	 marked_bams.view { "Mark Duplicates Output: $it" }

    // **Step 6: Split N CIGAR Reads**
    split_bams = SPLIT_NCIGAR_READS(marked_bams.map { 
            tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index, genome_dict)

    // **Step 7: SAMTOOLS CALMD**
    calmd_bams = SAMTOOLS_CALMD(split_bams.map { tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index)

    // **Step 8: Base Quality Score Recalibration (BQSR)**
    recalibrated_bams = GATK_RECALIBRATION(calmd_bams.map { 
            tuple(it[0], it[1], it[2], it[3]) }, reference_genome, genome_index, genome_dict, merged_vcf, merged_vcf_index)

    // **Step 9: Convert BED to Interval List**
    interval_list_ch = BED_TO_INTERVAL_LIST(denylist_bed, reference_genome, genome_dict)

    // **Step 10: Scatter the Interval List**
    scattered_intervals_ch = SCATTER_INTERVAL_LIST(interval_list_ch, genome_dict)

    // **Step 11: Perform Variant Calling using GATK HaplotypeCaller**
    gvcf_output = GATK_HAPLOTYPE_CALLER(recalibrated_bams.map {
            tuple(it[0], it[1], it[2], it[3]) 
        }, reference_genome, genome_index, genome_dict, scattered_intervals_ch)
	
	gvcf_output.view { "Raw GVCF output: $it" }
	
	// **Step 12: Generate variant statistics**
    bcftools_stats_ch = BCFTOOLS_STATS(gvcf_output)
	
    // **Step 13: Apply GATK Variant Filtering**
    filtered_individual_vcfs = GATK_VARIANT_FILTER(gvcf_output, reference_genome, genome_index, genome_dict)

    

    // **Step 14: Handle Individual vs Merged VCFs**
    if (params.merge_vcf) {
        log.info "🔄 Merging VCF files..."
        merged_filtered_vcf = BCFTOOLS_MERGE(filtered_individual_vcfs.map { it[1] }.collect())
        final_vcf_output = merged_filtered_vcf
    } else {
        log.info "📄 Keeping individual VCFs..."
        final_vcf_output = filtered_individual_vcfs
    }

    log.info "✅ Variant Calling Completed."

    emit:
        final_vcf = final_vcf_output
        stats_report = bcftools_stats_ch
}
