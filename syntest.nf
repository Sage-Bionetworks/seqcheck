process get_fastqs {
    output:
    stdout query_result

    """
    synapse -s query "select id from file where parentId=='syn8223228' and fileFormat=='fastq'" \
        | grep syn
    """

}

query_result
    .splitText()
    .map { it.replaceFirst(/\n/, '') }
    .into { fastq_entities }

process get_name {
    tag "$syn_id"

    input:
    val syn_id from fastq_entities

    output:
    stdout into syn_names

    """
    synapse -s query "select name from file where id=='${syn_id}'" \
        | awk '{if (NR>1) {print \$1}}'
    """
}

syn_names.subscribe { print it }
