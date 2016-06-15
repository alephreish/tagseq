Tag-Seq tools
=============

Included in this project are tools for Tag-Seq index design and data processing. The scripts `demultiplex_tags.pl` and `deduplicate_tags.pl` work in couple to allow demultiplexing and deduplication of Tag-Seq reads labelled with (moderately) degenerate barcode regions (mDBRs). See below on how to parallelize the processing of the data.

demultiplex_tags.pl
-------------------

A tool for demultiplexing Tag-Seq data from a sequencing run.

Use as: `demultiplex_tags.pl [arguments] <reads.fq(.gz)> <specimens.tsv>`

Input files:

* `reads.fq(.gz)`: fastq file (optionally gzipped)
* `specimens.tsv`: a list of specimen names with associated P5 and P7 index sequences
- it is expected to be a tab-separated file with at least three columns: name, i7\_seq and i5\_seq (in any order, headers must be specified on the first line)

Arguments:

* `-i`: Index read length [defaults to 8]
* `-y`: allowed mismatches when matching P7 index [1]
* `-o`: Output folder [.]
* `-Z`: don't gZip output files (switch)
* `-d`: Dry run (no output) (switch)

deduplicate_tags.pl
-------------------

A tool for deduplicating Tag-Seq reads from an individual demultiplexed library.

Use as: `deduplicate_tags.pl [arguments] <reads.fq(.gz)>`

Input files:

* `reads.fq(.gz)`: fastq file corresponding to data from one sample (optionally gzipped)

Arguments:

* `m`: Matching length [defaults to 10]
* `i`: Index read length [8]
* `x`: allowed mismatches when comparing mDBRs for PCR-duplication detection [1]
* `t`: Trim this number of bases from the beginning of every sequence [2]
* `o`: Output file [- (STDOUT)]
* `f`: use only this Fraction of the reads [1]
* `d`: Dry run (no output) (switch)

Parallelizing read processing
-----------------------------

The first step, demultiplexing, is parallelizable only if more than one fastq file is to be demultiplexed. The second step, deduplication, is easily run as multiple processes or jobs with each job working on one individual library. E.g.:

    ncpus=18
    mkdir -p demult dedup
    demultiplex_tags.pl -o demult reads.fq.gz specimens.tsv
    for demult in demult/*.fq.gz; do
        while [ $(jobs | wc -l) -ge "$ncpus" ]; do sleep 1; done
        lib=$(basename "$demult")
        deduplicate_tags.pl -t2 -o "dedup/$lib" "$demult" &
    done
    wait

mdbr.pl
-------

A script to generate groups of compatible mDBRs of a given length.

Depends on [cliquer](http://users.aalto.fi/~pat/cliquer.html) [available from Ubuntu repositories].

Use as: `perl mdbr.pl <DBR length> <difference threshold> <alphabet>`

Arguments:

* `DBR length`: length of the coding part of the indices
* `difference threshold`: min number of differing positions required for a pair of mDBRs to be considered as compatible
* `alphabet`: currently supported alphabets are `RY` and `ATGC`

Example:

    $ perl mdbr.pl 6 3 RY
    Searching for all maximum weight cliques...
    RRYRYY RRYYRR RYRRRR RYRYYY YRRRRY YRRYYR YYYRYR YYYYRY
    RRYRYY RRYYRR RYRRRR RYRYYY YRRRYR YRRYRY YYYRRY YYYYYR
    ....
    Clique size: 8; number of cliques: 240

polya\_trim.pl
--------------

A simple script to trim reads to get rid of: polyA tails, occasional polyG artifacts (NextSeq-specific), 5' adaptor overhangs, reads with low-quality indices. For advanced trimming options it is recommended to use [cutadapt](http://cutadapt.readthedocs.org/).

Use as:

    polya_trim.pl [arguments] < input.fq > output.fq

or

    gzip -cd input.fq.gz | polya_trim.pl [arguments] | gzip > output.fq.fq


Arguments:

* `L`: minimum read Length before or after trimming [default: 77]
* `A`: trim polyA-tails of at least this length     [1]
* `G`: trim polyG-tails of at least this length     [1]
* `5`: trim this Number of bases from 5' end* ` [0]
* `I`: discard reads with polyN and 'slipped' indexes (flag)

Example:

	$ gzip -cd myReads.fq.gz | perl polya_trim.pl -L 75 -A 10 -G 0 -5 2 -I | gzip -c > myTrimmedReads.fq.gz

Citation
--------

* Rozenberg A, Leese F, Weiss LC, Tollrian R (2016). *BioTechniques* 61, in print
