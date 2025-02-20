nextflow.enable.dsl = 2

// ✅ Include required annotation processes
include { ANNOTATE_INDIVIDUAL_VARIANTS } from '../modules/snpeff_annotate.nf'
include { ANNOTATE_INDIVIDUAL_VARIANTS_VEP } from '../modules/ensemblvep_annotate.nf'
include { ANNOTATE_VARIANTS } from '../modules/snpeff_annotate.nf'
include { ANNOTATEVARIANTS_VEP } from '../modules/ensemblvep_annotate.nf'
include { EXTRACT_VCF } from '../modules/extract_vcf.nf'
include { EXTRACT_individual_VCF } from '../modules/extract_individual_vcf.nf'

workflow ANNOTATION {

    take:
        final_vcf          // ✅ Input: Either individual or merged VCF output from VARIANT_CALLING
        merge_vcf          // ✅ Input: Flag indicating whether merged VCFs should be used
        snpeff_jar
        snpeff_config
        snpeff_db
        genomedb
        vep_cache
        clinvar_vcf
        clinvar_index

    main:
        log.info "🔬 Starting Variant Annotation Workflow..."

        if (merge_vcf) {
            log.info "🧬 Annotating Merged VCF..."

            // **Step 1: Annotate merged VCF with SnpEff**
            annotated_merged_vcf = ANNOTATE_VARIANTS(
                final_vcf,
                snpeff_jar,
                snpeff_config,
                snpeff_db,
                genomedb
            )

            // **Step 2: Annotate merged VCF with Ensembl VEP**
            annotated_merged_vcf_vep = ANNOTATEVARIANTS_VEP(
                final_vcf,  
                vep_cache,
                clinvar_vcf,
                clinvar_index
            )

            // **Step 3: Extract data from merged VCF**
            extracted_csv = EXTRACT_VCF(annotated_merged_vcf)

            log.info "✅ Merged VCF Annotation Completed."

        } else {
            log.info "🧬 Annotating Individual VCFs..."

            // **Step 1: Annotate individual VCFs with SnpEff**
            annotated_individual_vcfs = ANNOTATE_INDIVIDUAL_VARIANTS(
                final_vcf,
                snpeff_jar,
                snpeff_config,
                snpeff_db,
                genomedb
            )

            // **Step 2: Annotate individual VCFs with Ensembl VEP**
            annotated_individual_vcf_vep = ANNOTATE_INDIVIDUAL_VARIANTS_VEP(
                final_vcf,  
                vep_cache,
                clinvar_vcf,
                clinvar_index
            )

            // **Step 3: Extract data from individual VCFs**
            extracted_csv = EXTRACT_individual_VCF(annotated_individual_vcfs)

            log.info "✅ Individual VCF Annotation Completed."
        }

        log.info "✅ Variant Annotation Workflow Completed."

    // ✅ Correct Placement of `emit` (AFTER the conditional blocks)
    emit:
        annotated_vcf = merge_vcf ? annotated_merged_vcf_vep : annotated_individual_vcf_vep
        annotation_reports = extracted_csv
}
