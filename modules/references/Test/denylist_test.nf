// ========================== Download Denylist File ========================== //
process CHECK_OR_DOWNLOAD_DENYLIST {
    tag "Check or Download Denylist BED File"
    container null
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    output:
    path "denylist.bed", emit: denylist

    script:
    """
    wget -q -O denylist.bed.gz ${params.denylist_download_url}

    # Check if the file is gzipped and unzip if necessary
    if file denylist.bed.gz | grep -q 'gzip'; then
        gunzip denylist.bed.gz
    else
        mv denylist.bed.gz denylist.bed
    fi
    """
}