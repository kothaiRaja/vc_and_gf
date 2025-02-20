process STAR_ALIGNMENT {
    tag { sample_id }

    container "https://depot.galaxyproject.org/singularity/star%3A2.7.11a--h0033a41_0"
    publishDir "${params.outdir}/multiqc_input", mode: "copy"

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