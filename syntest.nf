process get_fastqs {
    output:
    stdout query_result

    """
    synapse -s query "select id,name from file where parentId=='syn8223228' and fileFormat=='fastq'"
    """

}

query_result
    .splitCsv(sep: '\t', header: true)
    // .subscribe { row ->
    //     println "${row['file.id']} - ${row['file.name']}"
    // }
    .into { fastq_entities }

process get_name {
    tag "${fastq_entity['file.id']}"

    input:
    val fastq_entity from fastq_entities

    output:
    stdout into syn_names
    val read_file into read_files

    script:
    syn_id = fastq_entity['file.id']
    read_file = "foo/${fastq_entity['file.name']}"
    """
    synapse -s query "select name from file where id=='${syn_id}'" \
        | awk '{if (NR>1) {print \$1}}'
    """
}

// syn_names.subscribe { print it }
read_files.subscribe { println it }
