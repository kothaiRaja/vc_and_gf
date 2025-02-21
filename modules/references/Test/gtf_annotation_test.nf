// ========================== Download GTF Annotation File ========================== //
process CHECK_OR_DOWNLOAD_GTF {
    tag "Check or Download GTF File"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output: 
    path "annotations.gtf", emit: gtf

    script:
    """
    wget -q -O annotations.gtf.gz ${params.gtf_download_url}

    # Check if the file is gzipped and unzip if necessary
    if file annotations.gtf.gz | grep -q 'gzip'; then
        gunzip annotations.gtf.gz
    else
        mv annotations.gtf.gz annotations.gtf
    fi
    """
}