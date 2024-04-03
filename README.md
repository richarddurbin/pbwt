The pbwt package provides a core implementation and development
environment for PBWT (Positional Burrows-Wheeler Transform) methods
for storing and computing on genome variation data sets.  

More precisely, PBWT supports a run-length compressed representation
of aligned haplotype data, on which efficient matching algorithms can
be built. Typically PBWT compression is much better than generic
compression, particularly for large numbers of haplotypes, and search
algorithms are linear in the query size independent of reference size.

A description of the basic data structure and matching algorithms is
given in "Efficient haplotype matching and storage using the
Positional Burrows-Wheeler Transform (PBWT)", Richard Durbin
[Bioinformatics 30:1266-72 (2014)](https://academic.oup.com/bioinformatics/article/30/9/1266/236397).

There are various phasing and imputation methods in the software that
are not yet published.

Richard Durbin <rd@sanger.ac.uk>

May 2013, updated September 2014

April 2024: the `-paint` and `-paintSparse` commands are described in Yaoling Yang, Richard Durbin, Astrid K. N. Iversen, and Daniel J. Lawson. 2024. “Sparse Haplotype-Based Fine-Scale Local Ancestry Inference at Scale Reveals Recent Selection on Immune Responses.” [medRxiv 2024.03.13.24304206](https://doi.org/10.1101/2024.03.13.24304206)


Installation instructions
------------------------
Download htslib from https://github.com/samtools/htslib, and compile it

   git clone https://github.com/samtools/htslib
   cd htslib
   make
   cd ..

Download and make pbwt

   git clone https://github.com/richarddurbin/pbwt
   cd pbwt
   make

Brief usage instructions
------------------------

Typing

    pbwt

by itself gives a list of commands with brief descriptions.

A quick synopsis for usage is:

    macs 11000 1e6 -t 0.001 -r 0.001 > 11k.macs
    pbwt -checkpoint 10000 -readMacs 11k.macs -write macs11k.pbwt -writeSites macs.sites

NB "checkpoint 10000" writes out files every 10000 sites during the vcfq 
conversion to alternating checkA.{pbwt,sites} and checkB.{pbwt,sites} files.

    pbwt -read macs11k.pbwt -subsample 0 10000 -write macs10k.pbwt
    pbwt -read macs11k.pbwt -subsample 10000 1000 -write macs1k.pbwt
    pbwt -read macs10k.pbwt -sfs > macs10k.sfs

gives the site frequency spectrum for macs10k

    pbwt -read macs1k.pbwt -haps macs1k.haps

writes out the haplotypes stored in macs1k

    pbwt -read macs10k.pbwt -matchDynamic macs1k.pbwt > macs1k-10k.max

for each sequence in macs1k, finds maximal matches to anything in macs10k

    pbwt -read macs10k.pbwt -maxWithin > macs10k.max

finds maximal matches for each sequence in macs10k to anything else in macs10k

To start from real data in a .vcf file rather than a macs simulation use

    pbwt -checkpoint 10000 -readVcfGT data.vcf -writeAll data

Note that -writeAll xxx will write xxx.pbwt, xxx.sites, xxx.samples and any other
associated files, and -readAll xxx will correspondingly read xxx.pbwt
and any available files based on suffix.

pbwt is very happy to handle up to 100,000 haplotypes, probably a
million.
