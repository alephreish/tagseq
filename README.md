Tag-Seq tools
=============

Included in this project are tools for Tag-Seq index design and data processing.

process_tags.pl
---------------

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

    > perl mdbr.pl 6 3 RY
    Searching for all maximum weight cliques...
    RRYRYY RRYYRR RYRRRR RYRYYY YRRRRY YRRYYR YYYRYR YYYYRY
    RRYRYY RRYYRR RYRRRR RYRYYY YRRRYR YRRYRY YYYRRY YYYYYR
    ....
    Clique size: 8; number of cliques: 240

polya\_trim.pl
--------------

A script to trim reads to get rid of: polyA tails, occasional polyG artifacts (NextSeq-specific), 5' adaptor overhangs, reads with low-quality indices

Use as: `perl polya_trim.pl [min length] [trim polyA] [trim polyG] [trim N bases from 5'] [trim indices] < [input fastq] > [output fastq]`

Arguments:

* `min length` (integer number): minimum read length before or after trimming, any read shorter than this number is discarded
* `trim polyA` (0 or 1): flag indicating whether to trim polyA tails
* `trim polyG` (0 or 1): flag indicating whether to trim ``polyG'' tails
* `trim N bases from 5'` (integer number): trim this number of bases from the 5'-end
* `trim indices` (0 or 1): flag, whether to remove reads with indices composed of polyN or the `AGATCTCG` ``slipping'' sequence

Example:

	>gzip -cd myReads.fq.gz | perl polya_trim.pl 75 1 0 2 1 | gzip -c > myTrimmedReads.fq.gz

Citation
--------

* Rozenberg A, Leese F, Weiss LC, Tollrian R (2016). *Biotechniques*, accepted
