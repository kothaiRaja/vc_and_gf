process DOWNLOAD_SNPEFF_TOOL {
    tag "Download SnpEff Tool"
    publishDir "${params.test_data_dir}/Tools", mode: 'copy'
    container null

    output:
    path "${params.snpeff_jar_dir}/snpEff.jar", emit: snpeff_jar
    path "${params.snpeff_jar_dir}/snpEff.config", emit: snpeff_config

    script:
    """
    mkdir -p ${params.snpeff_jar_dir}
    wget -q -O snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
    unzip -j snpEff_latest_core.zip -d ${params.snpeff_jar_dir}
    rm snpEff_latest_core.zip
    """
}

process DOWNLOAD_SNPEFF_DB {
    tag "Download SnpEff Database"
	publishDir "${params.test_data_dir}/Tools/snpEff", mode: 'copy'
	container null
	
	
    input:
    val genome
    path snpeff_jar_path

    output:
    path "${params.snpeff_db_dir}/${genome}"
	
	script:
    """
    
    mkdir -p ${params.snpeff_db_dir}

    
    data_dir=\$(realpath ${params.snpeff_db_dir})

    # Download the database
    # Check if the database is already downloaded
    if [ ! -d "\$data_dir/${genome}" ]; then
        echo "Downloading SnpEff database for ${genome}..."
        java -Xmx4g -Xms2g -jar ${snpeff_jar_path} download ${genome} -dataDir \$data_dir -v
    else
        echo "SnpEff database for ${genome} already exists. Skipping download."
    fi
    """
}

process DOWNLOAD_ARRIBA {
    container null
	publishDir "${params.test_data_dir}/Tools/ARRIBA", mode: 'copy' 
    
	output:
    path 'arriba_v2.4.0', emit: 'arriba_dir'
	
	script:
    """
   
    URL="https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz"
    TARGET_DIR="arriba_v2.4.0"

    # Download Arriba tarball
    wget -O arriba.tar.gz \$URL

    # Create the target directory and extract the tarball
    mkdir -p \$TARGET_DIR
    tar -xzvf arriba.tar.gz -C \$TARGET_DIR --strip-components=1

    # Clean up the downloaded tarball
    rm arriba.tar.gz

    # Output the extracted directory
    echo "Arriba tool extracted to: \$TARGET_DIR"
    """
}

