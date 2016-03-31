Tag-Seq tools
=============

Included in this project are tools for Tag-Seq index design and data processing.

process_tags.pl
---------------

A tools for demultiplexing Tag-Seq reads tagged with mDBRs.

Use as: `process_tags.pl [arguments] <reads.fq(.gz)> <specimens.txt>`

Input files:

* `reads.fq(.gz)`: the fastq file with the sequences (optionally gzipped)
* `specimens.txt`: a list of specimen names with associated P5 and P7 index sequences

Arguments:

* `-m`: Matching length [defaults to 10]
* `-i`: Index read length [8]
* `-x`: allowed mismatches when comparing mDBRs for PCR-duplication detection [1]
* `-y`: allowed mismatches when matching P7 index [1]
* `-t`: Trim this number of bases from the beginning of every sequence [2]
* `-z`: gZip output files (y/n) [y]
* `-d`: Dry run (no output) (y/n) [n]
* `-o`: Output folder [.]
* `-f`: use only this Fraction of the reads [all]


mdbr.pl
-------

A script to generate groups of compatible mDBRs of a given length.

Depends on [cliquer](http://users.aalto.fi/~pat/cliquer.html) [available from Ubuntu repositories].

Use as: `perl mdbr.pl <DBR length> <difference threshold> <alphabet>`

Arguments:

* `DBR length`: length of the coding part of the indices
* `difference threshold`: min number of differing positions required for a pair of mDBRs to be considered as compatible
* `alphabet`: currently supported alphabets are `RY` or `ATGC`

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
* `5`: trim this Number of bases from 5' end        [0]
* `I`: discard reads with polyN and 'slipped' indexes (flag)

Example:

	$ gzip -cd myReads.fq.gz | perl polya_trim.pl -L 75 -A 10 -G 0 -5 2 -I | gzip -c > myTrimmedReads.fq.gz

Citation
--------

* Rozenberg A, Leese F, Weiss LC, Tollrian R (2016). *Biotechniques*, accepted
