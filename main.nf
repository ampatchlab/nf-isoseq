#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
 *
 * ampatchlab/nf-isoseq: Nextflow pipeline for IsoSeq v3.1
 *
 * Copyright (C) 2019 QIMR Berghofer Medical Research Institute
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


nextflow.preview.dsl=2


import nextflow.splitter.CsvSplitter
import nextflow.config.ConfigParser

nextflow_config = file("${baseDir}/nextflow.config").text
parsed_config = new ConfigParser().setIgnoreIncludes(true).parse(nextflow_config)
defaults = parsed_config.params

check_params()



/*
 * Log workflow options
 */

// Nextflow exectution options
log.info("Config profile: ${workflow.profile}")
log.info("Workflow revision: ${workflow.revision}")
log.info("Workflow work-dir: ${workflow.workDir}")

// Input options
log.info("Input csv: ${params.csv}")

// Reference genome options
log.info("Reference genome: ${params.genome}")
log.info("Reference FASTA file: ${params.fasta}")
log.info("Reference GTF file: ${params.gtf}")

// Primers
log.info("Primers to trim: ${params.primers}")
log.info("Primers FASTA file: ${params.primer_seqs}")

// Output options
log.info("Output directory: ${params.outdir}")



/*
 * Log advanced options
 */

// Subread filtering
log.info("Min subread quality: ${params.min_subread_quality}")
log.info("Min subread length: ${params.min_subread_length}")
log.info("Skip subread filtering: ${String.valueOf(params.skip_subread_filtering)}")

// CCS calling params
log.info("Min passes for CCS calling: ${params.min_ccs_passes}")

// IsoSeq refine params
log.info("Require poly-A tails for refinement: ${String.valueOf(params.require_polya)}")

// IsoSeq cluster params
log.info("Max CCS reads for clustering: ${params.max_ccs_reads}")

// IsoSeq polish params
log.info("Read quality cutoff for polishing: ${params.rq_cutoff}")

// Reports options
log.info("Execution report: ${params.execution_report}")
log.info("Trace report: ${params.trace_report}")
log.info("Timeline report: ${params.timeline_report}")
log.info("Flowchart: ${params.flowchart}")

// AWS Batch options
log.info("AWS Batch JobQueue: ${params.aws_queue}")
log.info("AWS Region: ${params.aws_region}")



/*
 * Validate libraries and basecall inputs
 */

validate_input_csv()



/*
 * Processes
 */


// STEP 1 - BAX to BAM
// Create subreads from the raw basecalls

process bax2bam {

    tag { movie }

    label 'bax2bam'

    input:
    tuple movie, path(bax1), path(bax2), path(bax3)

    output:
    tuple movie, path("${movie}.subreads.bam{,.pbi}")

    """
    bax2bam -o "${movie}" "${bax1}" "${bax2}" "${bax3}"
    """
}


// STEP 2 - Filter subreads
// Filter the input subreads by read-quality and read-length

process create_filtered_subreadsets {

    tag { movie }

    label 'pbcoretools'

    input:
    tuple movie, path(indexed_subreads)

    output:
    tuple movie, path("${movie}.filtered.subreadset.xml")

    script:
    def (bam, pbi) = indexed_subreads

    """
    dataset create \\
        --type SubreadSet \\
        --relative \\
        "${movie}.subreadset.xml" \\
        "${bam}"
    dataset filter \\
        "${movie}.subreadset.xml" \\
        "${movie}.filtered.subreadset.xml" \\
        'rq>=${params.min_subread_quality}' \\
        'length>=${params.min_subread_length}'
    dataset relativize "${movie}.filtered.subreadset.xml"
    """
}

process apply_subread_filters {

    tag { movie }

    label 'pbbam'

    input:
    tuple movie, path(subreadset), path(indexed_subreads)

    output:
    tuple movie, path("${movie}.filtered.subreads.bam{,.pbi}")

    """
    pbmerge \\
        -o "${movie}.filtered.subreads.bam" \\
        "${subreadset}"
    """
}


// STEP 3 - Subreads to CCS
// Call circular consensus sequences

process ccs {

    tag { movie }

    label 'pbccs'

    input:
    tuple movie, path(indexed_subreads)

    output:
    tuple movie, path("${movie}.ccs.bam{,.pbi}"), emit: bam
    path "${movie}.ccs_report.txt", emit: report

    script:
    def (bam, pbi) = indexed_subreads

    """
    ccs \\
        --numThreads ${task.cpus} \\
        --minPasses "${params.min_ccs_passes}" \\
        --noPolish \\
        --reportFile "${movie}.ccs_report.txt" \\
        "${bam}" \\
        "${movie}.ccs.bam"
    """
}


// STEP 4 - CCS to FL
// Generate full-length reads by primer removal and demultiplexing

process lima {

    tag { movie }

    label 'lima'

    input:
    tuple movie, path(indexed_ccs)
    tuple val(five_primer), val(three_primer)
    path barcodes

    output:
    tuple movie, path("${movie}.fl.${five_primer}--${three_primer}.bam{,.pbi}"), emit: bam
    path "${movie}.fl.lima.counts", emit: counts
    path "${movie}.fl.lima.report", emit: report
    path "${movie}.fl.lima.summary", emit: summary

    script:
    def (bam, pbi) = indexed_ccs

    """
    lima \\
        --isoseq \\
        --num-threads ${task.cpus} \\
        "${bam}" \\
        "${barcodes}" \\
        "${movie}.fl.bam"
    """
}


// STEP 5 - FL to FLNC
// Remove concatemers and optionally poly-A tails

process refine {

    tag { movie }

    label 'isoseq3'

    input:
    tuple movie, path(indexed_ccs)
    path barcodes

    output:
    tuple movie, path("${movie}.flnc.bam{,.pbi}"), emit: bam
    path "${movie}.flnc.report.csv", emit: report

    script:
    def require_polya = params.require_polya ? '--require-polya' : ''
    def (bam, pbi) = indexed_ccs

    """
    isoseq3 refine \\
        --verbose \\
        ${require_polya} \\
        "${bam}" \\
        "${barcodes}" \\
        "${movie}.flnc.bam"
    """
}


// STEP 6 - FLNC to UNPOLISHED
// Cluster FLNC reads and generate unpolished transcripts

process create_transcriptsets {

    tag { library }

    label 'pbcoretools'

    input:
    tuple library, path(bams), path(indexes)

    output:
    tuple library, path("${library}.transcriptset.xml")

    script:
    def quoted_bams = bams.collect { /"${it}"/ }.join(' ')

    """
    dataset create \\
        --type TranscriptSet \\
        --relative \\
        "${library}.transcriptset.xml" \\
        ${quoted_bams}
    """
}

process cluster {

    tag { library }

    label 'isoseq3'

    input:
    tuple library, path(transcriptset), path(bams), path(indexes)

    output:
    tuple library, path("${library}.unpolished.bam{,.pbi}"), emit: bam

    """
    isoseq3 cluster \\
        --num-threads ${task.cpus} \\
        --verbose \\
        --poa-cov "${params.max_ccs_reads}" \\
        "${transcriptset}" \\
        "${library}.unpolished.bam"
    """
}


// STEP 7 - UNPOLISHED to POLISHED
// Polish transcripts using subreads

process create_subreadsets {

    tag { library }

    label 'pbcoretools'

    input:
    tuple library, path(bams), path(indexes)

    output:
    tuple library, path("${library}.subreadset.xml")

    script:
    def quoted_bams = bams.collect { /"${it}"/ }.join(' ')

    """
    dataset create \\
        --type SubreadSet \\
        --relative \\
        "${library}.subreadset.xml" \\
        ${quoted_bams}
    """
}

process polish {

    tag { library }

    label 'isoseq3'

    input:
    val library
    tuple path(unpolished), path(unpolished_index)
    tuple path(subreadset), path(subreads), path(subreads_indexes)

    output:
    tuple library, path("${library}.polished.bam{,.pbi}"), emit: bam
    tuple library, path("${library}.polished.hq.fasta.gz"), emit: hq_reads
    tuple library, path("${library}.polished.lq.fasta.gz"), emit: lq_reads
    path "${library}.polished.cluster_report.csv", emit: report

    """
    isoseq3 polish \\
        --num-threads ${task.cpus} \\
        --verbose \\
        --rq-cutoff "${params.rq_cutoff}" \\
        --coverage "${params.max_subreads}" \\
        "${unpolished}" \\
        "${subreadset}" \\
        "${library}.polished.bam"
    """
}


// STEP 8 - POLISHED to SAM
// Align the polished transcripts using Minimap2

process minimap2 {

    tag { library }

    label 'minimap2'

    input:
    tuple library, path(reads)
    path ref

    output:
    file "${library}.sam"

    """
    minimap2 \\
        -t ${task.cpus} \\
        -ax splice:hq \\
        -uf \\
        -o "${library}.sam" \\
        "${ref}" \\
        "${reads}"
    """
}


// STEP 9 - SAM to BAM
// Create BAM and index files for IGV review

process samtools {

    tag { sam.getBaseName() }

    label 'samtools'

    input:
    path sam

    output:
    path "${sam.getBaseName()}.bam"

    """
    samtools view -u "${sam}" | samtools sort -o "${sam.getBaseName()}.bam" -
    samtools index "${sam.getBaseName()}.bam"
    """
}


// STEP 10 - SAM to Junctions TSV
// Evaluate splice junction consistency with known annotations

process junceval {

    tag { sam.getBaseName() }

    label 'minimap2'

    input:
    path sam
    path gtf

    output:
    path "${sam.getBaseName()}.junctions.tsv"

    """
    paftools.js junceval \\
        -p \\
        "${gtf}" \\
        "${sam}" \\
        > "${sam.getBaseName()}.junctions.tsv"
    """
}



/*
 * Functions
 */

def parse_input_csv( csv ) {

    def rows = Channel.fromPath( csv ) | splitCsv( header:true )

    def inputs = [:]

    inputs.libraries = rows \
        | map { row ->
            def movie = row.movie ?: file(row.bax1).getSimpleName()
            tuple(row.library.replaceAll(/\./, '_'), movie)
        }

    inputs.basecalls = rows \
        | map { row ->
            def movie = row.movie ?: file(row.bax1).getSimpleName()
            tuple(movie, file(row.bax1), file(row.bax2), file(row.bax3))
        }

    return inputs
}

def get_primer_labels( primers ) {

    labels = Channel.fromPath( primers )
        .splitFasta( record: [id: true, seqString: false ], limit: 2 )
        .collect { record -> record.id }

    return labels
}

def group_reads_by_library( reads, libraries ) {

    grouped = libraries \
        | groupTuple() \
        | map { library, movies -> tuple( groupKey(library, movies.size()), movies ) } \
        | transpose() \
        | map { library, movie -> tuple(movie, library) } \
        | join( reads ) \
        | map { movie, library, bam -> tuple(library, *bam) } \
        | groupTuple() \
        | map { library, bams, indexes -> tuple(library.toString(), bams, indexes) }

    return grouped
}

def usage() {

    log.info"""
    Usage:
        nextflow run -profile <profile> -revision <revision> ampatchlab/nf-isoseq [options]


    Nextflow execution options:

        -profile STR
            Nextflow configuration profile to use. Available profiles include:
            'awsbatch', 'conda', 'docker' and 'singularity'

        -revision STR
            Git branch/tag (version) of this workflow to use

        -work-dir DIR
            Directory where intermediate result files are stored

        -help
            Show additional execution options and exit


    Input options:

        --csv FILE
            Comma-separated list of libraries and basecall inputs


    Reference genome options:

        --genome STR
            Reference genome name [Either: ${defaults.genomes.keySet().join(", ")}; Default: ${defaults.genome}]

        --fasta FILE
            Override the reference genome FASTA with FILE [Default: ${defaults.fasta ?: null}]

        --gtf FILE
            Override the reference genome GTF with FILE [Default: ${defaults.gtf ?: null}]


    Sequencing primer options:

        --primers STR
            The primers to trim [Either: ${defaults.primer_sets.keySet().join(", ")}; Default: ${defaults.primers}]

        --primer_seqs FILE
            Override the primer sequences with those in FILE [Default: ${defaults.primer_seqs ?: null}]


    Output options:

        --outdir DIR
            Path where the results will be saved [Default: ${defaults.outdir}]


    Standard options:

        --advanced
            Show advanced usage and exit

        --help
            Show this message and exit

        --version
            Show the pipeline version and exit
    """.stripIndent()
}

def advanced() {

    log.info"""
    Subread filtering options:

        --min_subread_quality FLOAT
            Minimum read score of input subreads [Default: ${defaults.min_subread_quality}]

        --min_subread_length INT
            Minimum read length of input subreads [Default: ${defaults.min_subread_length}]

        --skip_subread_filtering
            Do not apply any of the above subread filters [Default: ${String.valueOf(defaults.skip_subread_filtering)}]


    CCS calling options:

        --min_ccs_passes INT
            Minimum number of subreads required to generate CCS [Default: ${defaults.min_ccs_passes}]


    IsoSeq refine options:

        --require_polya
            Require FL reads to have a poly(A) tail and remove it [Default: ${String.valueOf(defaults.require_polya)}]


    IsoSeq cluster options:

        --max_ccs_reads INT
            Maximum number of CCS reads used for POA consensus [Default: ${defaults.max_ccs_reads}]


    IsoSeq polish options:

        --rq_cutoff FLOAT
            Read quality cutoff for fastx output [Default: ${defaults.rq_cutoff}]

        --max_subreads INT
            Maximum number of subreads used for polishing [Default: ${defaults.max_subreads}]


    Report options:

        --execution_report STR
            Name of the Nextflow execution report to generate [Default: ${defaults.execution_report}]

        --trace_report STR
            Name of the Nextflow trace report to generate [Default: ${defaults.trace_report}]

        --timeline_report STR
            Name of the Nextflow timeline report to generate [Default: ${defaults.timeline_report}]

        --flowchart STR
            Name of the Nextflow flowchart to generate [Default: ${defaults.flowchart}]


    AWS Batch options:

        --aws_queue STR
            AWS Batch JobQueue definition [Default: ${defaults.aws_queue}]

        --aws_region STR
            AWS Region definition [Default: ${defaults.aws_region}]
    """.stripIndent()
}

def die() {

    usage()
    exit 1
}

def check_params() {

    // Standard options

    if (params.advanced) {
        advanced()
        exit 0
    }

    if (params.help) {
        usage()
        exit 0
    }

    if (params.version) {
        log.info(workflow.manifest.version)
        exit 0
    }


    // Required options

    if (!params.csv) {
        log.error("A list of libraries and basecalls is required. Please use the `--csv` option.")
        die()
    }

    if (file(params.csv).getExtension() != "csv") {
        log.error("CSV input file `${params.csv}` must be a CSV file with the '.csv' extension.")
        die()
    }


    // Reference genome options

    if (!params.genome) {
        log.error("Please specify a value for `--genome`; can be one of ${params.genomes.keySet().join(", ")}")
        die()
    }

    params.fasta = params.genome ? params.genomes[ params.genome ].fasta : null
    params.gtf = params.genome ? params.genomes[ params.genome ].gtf : null

    if (!params.fasta) {
        log.error("A reference FASTA file is required. Please use the `--fasta` option.")
        die()
    }

    if (!params.gtf) {
        log.error("A reference GTF file is required. Please use the `--gtf` option.")
        die()
    }


    // Sequencing primers

    if (!params.primers) {
        log.error("Please specify a value for `--primers`; can be one of ${params.primer_sets.keySet().join(", ")}")
        die()
    }

    params.primer_seqs = params.primer_sets ? params.primer_sets[ params.primers ].fasta : null

    if (!params.primer_seqs) {
        log.error("A pair of primer sequences is required. Please use the `--primer_seqs` option.")
        die()
    }


    // Report options

    if (!params.execution_report.toString().endsWith('.html')) {
        log.error("The filename specified using `--execution_report` must end with '.html'")
        die()
    }

    if (!params.trace_report.toString().endsWith('.txt')) {
        log.error("The filename specified using `--trace_report` must end with '.txt'")
        die()
    }

    if (!params.timeline_report.toString().endsWith('.html')) {
        log.error("The filename specified using `--timeline_report` must end with '.html'")
        die()
    }

    def flowchart_extns = ['.dot', '.html', '.pdf', '.png', '.svg']

    if (!(flowchart_extns.any { params.flowchart.toString().endsWith(it) })) {
        log.error("The filename specified using `--flowchart` must end with one of ${flowchart_extns.join(", ")}")
        die()
    }
}

def validate_input_csv() {

    def csv = new CsvSplitter().target(file(params.csv)).options(header:true)
    def rows = csv.list()

    def basecall_columns = ["bax1", "bax2", "bax3"]
    def required_columns = ["library"] + basecall_columns

    required_columns.each { col ->

        if (!csv.columnsHeader.contains(col)) {
            log.error("CSV input file `${params.csv}` does not contain a '${col}' column. Exiting.")
            exit 1
        }
    }

    def basecalls = []

    log.info("Validating ${rows.size()} entries...")

    rows.indexed(1).each { idx, row ->

        log.info("Validating entry ${idx}...")

        if (!row.library) {
            log.error("Entry ${idx}: Invalid 'library' value. Exiting.")
            exit 1
        }

        if (basecall_columns.collect { file(row."${it}").getSimpleName() }.toSet().size() != 1) {
            log.error("Entry ${idx}: Could not extract a unique movie name from the supplied bax files. Exiting.")
            exit 1
        }

        for (col in basecall_columns.indexed(1)) {
            bax = file(row."${col.value}")

            if (! bax.getName().endsWith(".${col.key}.bax.h5")) {
                log.error("Entry ${idx}: Invalid bax file extension: ${bax.getName()}")
                exit 1
            }

            if (bax.getSimpleName() + ".${col.key}.bax.h5" != bax.getName()) {
                log.error("Entry ${idx}: Invalid bax filename: ${bax.getName()}")
                exit 1
            }

            if (bax.getName() in basecalls) {
                log.error("Entry ${idx}: File '${bax.getName()}' cannot be used more than once. Exiting.")
                exit 1
            }

            basecalls.add(bax.getName())
        }
    }

    log.info("Done")
}



/*
 * Workflows
 */

workflow {

    main:
    input = parse_input_csv(params.csv)

    // bax2bam
    bax2bam(input.basecalls)

    // subreads
    subreads = ! params.skip_subread_filtering
             ? create_filtered_subreadsets(bax2bam.out) | join(bax2bam.out) | apply_subread_filters
             : bax2bam.out

    // consensus
    ccs(subreads)

    // demultiplex
    primer_ids = get_primer_labels(params.primer_seqs)
    lima(ccs.out.bam, primer_ids, params.primer_seqs)

    // refine
    refine(lima.out.bam, params.primer_seqs)

    // cluster
    grouped_flnc_reads = group_reads_by_library(refine.out.bam, input.libraries)
    create_transcriptsets(grouped_flnc_reads) \
        | join( grouped_flnc_reads) \
        | cluster

    // polish
    grouped_subreads = group_reads_by_library(subreads, input.libraries)
    create_subreadsets(grouped_subreads) \
        | join( grouped_subreads ) \
        | join( cluster.out.bam ) \
        | fork { library, subreadset, subreads, indexes, transcripts ->
            libraries: library
            transcripts: transcripts
            subreads: tuple(subreadset, subreads, indexes)
        } \
        | polish

    // minimap2
    minimap2(polish.out.hq_reads, params.fasta) | samtools

    // paftools
    junceval(minimap2.out, params.gtf)

    publish:
    subreads to: "${params.outdir}/subreads", mode: 'copy'
    ccs.out to: "${params.outdir}/ccs", mode: 'copy'
    lima.out to: "${params.outdir}/lima", mode: 'copy'
    refine.out to: "${params.outdir}/refine", mode: 'copy'
    cluster.out to: "${params.outdir}/cluster", mode: 'copy'
    polish.out to: "${params.outdir}/polish", mode: 'copy'
    minimap2.out to: "${params.outdir}/minimap2", mode: 'copy'
    samtools.out to: "${params.outdir}/samtools", mode: 'copy'
    junceval.out to: "${params.outdir}/junceval", mode: 'copy'
}
