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
Positional Burrows-Wheeler Transform (PBWT)", Richard Durbin (2013),
in preparation.

Phasing and imputation methods are currently under development.

Richard Durbin <rd@sanger.ac.uk>, May 2013


Installation
------------

The required `HTSlib` is bundled as a git submodule. To obtain `pbwt`, use
`git clone --recursive https://github.com/richarddurbin/pbwt.git`. This will
obtain both `pbwt` and `HTSlib` for you. Then, the following commands will
compile both `HTSlib` and `pbwt` for you:

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

    pbwt -read macs1k -haps macs1k.haps

writes out the haplotypes stored in macs1k

    pbwt -read macs10k.pbwt -maxIndexed macs1k.pbwt > macs1k-10k.max

for each sequence in macs1k, finds maximal matches to anything in macs10k

    pbwt -read macs10k.pbwt -maxWithin > macs10k.max

finds maximal matches for each sequence in macs10k to anything else in macs10k

To start from real data in a .vcf file rather than a macs simulation
we recommend to use the htscmd program (available from https://github.com/samtools/htslib.git)
to extract a simplified .vcfq file, e.g.

    htscmd vcfquery -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' < data.vcf > data.vcfq
    pbwt -checkpoint 10000 -vcfq data.vcfq -write data.pbwt -writeSites data.sites

pbwt is very happy to handle up to 100,000 haplotypes, probably a
million.
