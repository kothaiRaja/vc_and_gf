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