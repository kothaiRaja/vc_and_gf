// ========================== Create Genome Dictionary ========================== //
process CREATE_GENOME_DICT {
    tag "Create Genome Dictionary"
    container "https://depot.galaxyproject.org/singularity/picard%3A2.27.4--hdfd78af_0"
    publishDir "${params.test_data_dir}/reference", mode: 'copy'

    input:
    path genome_fa

    output:
    path "genome.dict", emit: genome_dict

    script:
    """
    echo " Creating genome dictionary using Picard..."
    picard CreateSequenceDictionary R=$genome_fa O=genome.dict
    """
}