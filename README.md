# ampatchlab/nf-isoseq

[![Build Status](https://codebuild.ap-southeast-2.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoidjdFOFhrK2FubFQ2LzlCaFROQVZ2Z3Qvczd0TjJXbzNhOXlKYTZsYUFIQUJMYStjOG1LL1p1cDA1Z2d2SEFFeFJjelBHS0kvM0VzWHZQSzZVUVlqT25BPSIsIml2UGFyYW1ldGVyU3BlYyI6IjBadStpd1ErMzhmSDVGbnYiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=master)](https://ap-southeast-2.console.aws.amazon.com/codesuite/codebuild/projects/nf-isoseq/history)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

[IsoSeq v3.1](https://github.com/PacificBiosciences/IsoSeq/blob/master/README_v3.1.md) RSII Nextflow pipeline

## Usage

```
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
        Reference genome name [Either: GRCh38, GRCm38; Default: GRCh38]

    --fasta FILE
        Override the reference genome FASTA with FILE [Default: null]

    --gtf FILE
        Override the reference genome GTF with FILE [Default: null]


Sequencing primer options:

    --primers STR
        The primers to trim [Either: Clontech, NEB; Default: Clontech]

    --primer_seqs FILE
        Override the primer sequences with those in FILE [Default: null]


Output options:

    --outdir DIR
        Path where the results will be saved [Default: ./results]


Standard options:

    --advanced
        Show advanced usage and exit

    --help
        Show this message and exit

    --version
        Show the pipeline version and exit
```

## Advanced options

```
Subread filtering options:

    --min_subread_quality FLOAT
        Minimum read score of input subreads [Default: 0.75]

    --min_subread_length INT
        Minimum read length of input subreads [Default: 35]

    --skip_subread_filtering
        Do not apply any of the above subread filters [Default: false]


CCS calling options:

    --min_ccs_passes INT
        Minimum number of subreads required to generate CCS [Default: 1]


IsoSeq refine options:

    --require_polya
        Require FL reads to have a poly(A) tail and remove it [Default: true]


IsoSeq cluster options:

    --max_ccs_reads INT
        Maximum number of CCS reads used for POA consensus [Default: 10]


IsoSeq polish options:

    --rq_cutoff FLOAT
        Read quality cutoff for fastx output [Default: 0.99]

    --max_subreads INT
        Maximum number of subreads used for polishing [Default: 60]


Report options:

    --execution_report STR
        Name of the Nextflow execution report to generate [Default: ./reports/execution_report.html]

    --trace_report STR
        Name of the Nextflow trace report to generate [Default: ./reports/trace_report.txt]

    --timeline_report STR
        Name of the Nextflow timeline report to generate [Default: ./reports/timeline_report.html]

    --flowchart STR
        Name of the Nextflow flowchart to generate [Default: ./reports/flowchart.png]


AWS Batch options:

    --aws_queue STR
        AWS Batch JobQueue definition [Default: false]

    --aws_region STR
        AWS Region definition [Default: false]
```

## Inputs

The input CSV must have the following required columns:

 * library: Unique library name or ID (required)
 * movie: Unique movie name or ID (optional)
 * bax1: Absolute path to movie.1.bax.h5 (required)
 * bax2: Absolute path to movie.2.bax.h5 (required)
 * bax3: Absolute path to movie.3.bax.h5 (required)

A single library may have multiple movies. Movies for each library will be merged during clustering.
Polished transcripts will then be aligned using [Minimap2](https://github.com/lh3/minimap2) and the
splice junctions evaluated using [paftools.js](https://github.com/lh3/minimap2/tree/master/misc).
