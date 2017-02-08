
version = 0.1

params.indir = false
params.reads = "${baseDir}/data/reads/*.fastq*"
params.build = false
params.mapper = "salmon"
params.index = params.build ? params.genomes[ params.build ][ params.mapper ] ?: false : false
params.fasta = params.build ? params.genomes[ params.build ].fasta ?: false : false
params.fasta_url = params.build ? params.genomes[ params.build ].fasta_url ?: false : false
params.save_reference = true
params.refdir = "${baseDir}/data/reference"
params.outdir = "${baseDir}/results"


log.info "============================================"
log.info "seqcheck : Sequence data sanity check  v${version}"
log.info "============================================"
log.info "Reads          : ${params.reads}"
log.info "Build          : ${params.build}"
log.info "Current home   : $HOME"
log.info "Current user   : $USER"
log.info "Current path   : $PWD"
log.info "Script dir     : $baseDir"
log.info "Working dir    : $workDir"
log.info "Reference dir  : ${params.refdir}"
log.info "index: ${params.index}"
log.info "fasta: ${params.fasta}"
log.info "Output dir     : ${params.outdir}"
log.info "============================================\n"

// Validate inputs
if( !params.build ){
    exit 1, "No reference genome specified."
}
else {
    index_file = file(params.index)
    if( !index_file.exists() ){
        log.info "Quant index not found: ${params.index}"
        fasta_file = file(params.fasta)
        if( !fasta_file.exists() ){
            log.info "FASTA file not found: ${params.fasta}"
            if( !params.download_fasta){
                exit 1, "No URL to download FASTA provided."
            }
            else {
                log.info "...will download from ${params.fasta_url}."
                run_download_fasta = true
            }
        }
        else {
            log.info "...will build index from ${params.fasta}."
            run_download_fasta = false
            fasta = file(params.fasta)
        }
    }
    else {
        run_download_fasta = false
        fasta = false
        quant_index = Channel
            .fromPath(params.index)
            .toList()
    }
}

/*
 * Create a channel for input read files
 */
Channel
    .fromPath( params.reads )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
    .into { read_files_fastqc; read_files_quant }

/*
 * PREPROCESSING - Download FASTA
 */
if( run_download_fasta ){
    log.info "\nDownloading FASTA from Ensembl..."
    process download_fasta {
        tag "${params.fasta_url}"
        publishDir path: "${params.refdir}/sequence/${params.build}", saveAs: { params.save_reference ? it : null }, mode: 'copy'

        output:
        file "*.{fa,fasta}" into fasta

        script:
        """
        curl -O -L ${params.download_fasta}
        if [ -f *.tar.gz ]; then
            tar xzf *.tar.gz
        elif [ -f *.gz ]; then
            gzip -d *.gz
        fi
        """
    }
}

/*
 * PREPROCESSING - Build quant index
 */
if( fasta ){
    if( params.mapper == "salmon" ){
        log.info "\nBuilding Salmon index from transcriptome FASTA..."
        process make_salmon_index {
            validExitStatus 0,137
            tag fasta
            publishDir path: "${params.refdir}/indexes/salmon", saveAs: { params.save_reference ? it : null }, mode: 'copy'

            input:
            file fasta from fasta

            output:
            file "${params.build}" into quant_index

            script:
            """
            mkdir ${params.build}
            salmon index \
                -t $fasta \
                -i ${params.build}
            """
        }
    }
}

/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$reads"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads -o .
    """
}

/*
 * STEP 2 - quantify
 */
if( params.mapper == "salmon" ){
    process salmon_quant {
        tag "$reads"
        publishDir "${params.outdir}/quant", mode: 'copy'

        input:
        file reads from read_files_quant
        file index from quant_index.first()

        output:
        file "${prefix}_quant" into quant_results

        script:
        prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
        """
        salmon quant \
            -i ${index} -l U \
            -r ${reads} \
            -o ${prefix}_quant
        """
    }
}

/*
 * STEP 11 MultiQC
 */
// process multiqc {
//     publishDir "${params.outdir}/MultiQC", mode: 'copy'
//
//     input:
//     file ('fastqc/*') from fastqc_results.flatten().toList()
//     file ('salmon/*') from alignment_logs.flatten().toList()
//
//     output:
//     file "*multiqc_report.html"
//     file "*multiqc_data"
//
//     script:
//     """
//     multiqc -f .
//     """
// }

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}\n"
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
