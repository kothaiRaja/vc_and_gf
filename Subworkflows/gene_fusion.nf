nextflow.enable.dsl = 2

// ✅ Include required processes
include { STAR_ALIGN_FUSION } from '../modules/star_align_fusion.nf'
include { ARRIBA } from '../modules/arriba.nf'
include { ARRIBA_VISUALIZATION } from '../modules/arriba_visualization.nf'

workflow GENE_FUSION {
    
    take:
        trimmed_reads_ch      // Trimmed reads from preprocessing
        star_genome_index     // STAR genome index
        reference_genome      // Reference genome FASTA
        gtf_annotation        // GTF annotation file
        arriba_blacklist      // ARRIBA blacklist file
        arriba_known_fusions  // ARRIBA known fusions file
        arriba_scripts_dir    // ARRIBA visualization scripts

    main:
        
        log.info "🧬 Starting Gene Fusion Detection Workflow..."

        // **Step 1: STAR Alignment for Fusion Detection**
        star_fusion_aligned = STAR_ALIGN_FUSION(trimmed_reads_ch, star_genome_index, gtf_annotation)

        // **Step 2: ARRIBA Fusion Detection**
        arriba_results = ARRIBA(
            star_fusion_aligned,
            reference_genome,
            gtf_annotation,
            arriba_blacklist,
            arriba_known_fusions
        )

        // **Step 3: ARRIBA Visualization**
        fusion_visuals = ARRIBA_VISUALIZATION(
            arriba_results,
            arriba_scripts_dir,
            reference_genome,
            gtf_annotation
        )

        log.info "✅ Gene Fusion Detection Workflow Completed."

    emit:
        fusion_results = arriba_results.fusions
		discarded_results = arriba_results.fusions_discarded
        fusion_visualizations = fusion_visuals
}
