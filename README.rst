dnapy
==========

.. image:: https://travis-ci.org/sherrillmix/dnapy.svg?branch=master
    :alt: Travis CI
    :target: https://travis-ci.org/sherrillmix/dnapy
.. image:: https://codecov.io/github/sherrillmix/dnapy/coverage.svg?branch=master
    :alt: Codecov
    :target: https://codecov.io/github/sherrillmix/dnapy?branch=master


Some python functions to deal with DNA alignments.
 
Installation
------------

Github
~~~~~~

To install the development version from github, clone the repository to a local directory using something like::

    git clone https://github.com/sherrillmix/dnapy.git

and run `setup.py` from the resulting directory (the `--user` installs it locally and doesn't require root access)::

  cd dnapy
  python setup.py install --user
  python setup.py test

To use the scripts directly from command line, e.g. `countbases file.bam`, (and to pass all tests above), you will need to make sure that your PATH contains the bin directory installed in by `setup.py` above. On Linux, this would mean putting something like::

   if [ -d "$HOME/.local/bin" ] ; then
       PATH="$HOME/.local/bin:$PATH"
   fi

in your `.bashrc` or `.profile`.

Usage
-----
The package provides the scripts:

abitotrace
~~~~

::
  
  usage: abitotrace [-h] abiFile
  
  A program to take an ABI .ab1 file from Sanger sequencing and output the
  traces to standard out as four columns.
  
  positional arguments:
    abiFile     an ABI .ab1 file containing the sequence read
  
  optional arguments:
    -h, --help  show this help message and exit
  
bamtoalign
~~~~

::
  
  usage: bamtoalign [-h] -s REFSEQ [-q MINQUALITY] [-v] [-r REGION] [-e ENDSPAN]
                    bamFile
  
  A program to convert a bam file into an aligned fasta file. The command
  generates fasta formatted output (two lines for each sequence: a name line
  prepended by > and a line containing the aligned sequence) to standard out.
  
  positional arguments:
    bamFile               a bam file containing the alignment
  
  optional arguments:
    -h, --help            show this help message and exit
    -s REFSEQ, --refseq REFSEQ
                          fasta file giving the reference sequence of interest
    -q MINQUALITY, --minQuality MINQUALITY
                          don't count alignments with a mapping quality less
                          than this
    -v, --verbose         increase output verbosity to stderr
    -r REGION, --region REGION
                          the region to pull reads from (note that the
                          underlying pysam does not like single base regions
                          like ch1:25. These instead be specified as
                          chr1:25-25.)
    -e ENDSPAN, --endSpan ENDSPAN
                          ignore spans of matches at the start or end of a read
                          less than this cutoff
  
countbases
~~~~

::
  
  usage: countbases [-h] [-v] [-r REGION] [-s] [-q MINQUALITY] bamFile
  
  A program to count the number of bases at each position in a region. The
  command generates standard output with columns referenceName, position,
  numberOfReads, and numbers of A, C, G, T (or A+, A-, C+, C-, G+, G-, T+, T- if
  --strand).
  
  positional arguments:
    bamFile               a bam file containing the alignment
  
  optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         increase output verbosity to stderr
    -r REGION, --region REGION
                          the region to count in
    -s, --strand          break base counts into positive and negative strand
                          alignments
    -q MINQUALITY, --minQuality MINQUALITY
                          don't count bases with a quality less than this
  
countkmers
~~~~

::
  
  usage: countkmers [-h] [-k KMERLENGTH] [-t NTHREADS]
                    fastqFiles [fastqFiles ...]
  
  A program to take a fastq file(s) and count the total k-mers across all reads
  in each file. Note that partial kmers are discarded e.g. the last 3 reads of a
  23 base read will be ignored. Return a comma separated file with a header row
  then a row for each file and a file column then a column for each kmer
  
  positional arguments:
    fastqFiles            a fastq file(s) (potentially gzipped) containing the
                          sequence reads
  
  optional arguments:
    -h, --help            show this help message and exit
    -k KMERLENGTH, --kmerLength KMERLENGTH
                          the lengh of kmer to be used. Be careful with values
                          larger than 20.
    -t NTHREADS, --nThreads NTHREADS
                          the number of threadss to use for processing. Should
                          be less than or equal the number of threads on
                          computer.
  
getstartends
~~~~

::
  
  usage: getstartends [-h] [-v] [-g MAXGAPS] [-r REGION] [-f FILE] [-n] [-c]
                      bamFile
  
  A program to pull start and end positions in a region. The command generates
  standard output with columns referenceName, start (1-based), end (1-based),
  strand
  
  positional arguments:
    bamFile               a bam file containing the alignment
  
  optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         increase output verbosity to stderr
    -g MAXGAPS, --maxGaps MAXGAPS
                          maximum allowed insertions or deletions in a read.
                          Otherwise discard
    -r REGION, --region REGION
                          the region to count in
    -f FILE, --file FILE  a text file specifying several regions to count where
                          each line gives a region e.g. chr1:1-100
    -n, --noHeader        suppress the initial header on csv output
    -c, --regionColumn    specify target region in first column (default: don't
                          show column)
  
removereads
~~~~

::
  
  usage: removereads [-h] [-d DOTS] [-o [OUTPUTFILES ...]] -f FILTERFILE [-k]
                     fastqFiles [fastqFiles ...]
  
  A program to filter reads by name from a single/set of fastq file(s). The
  script looks for reads which have a name line where the string before a space
  exactly matches a pattern. If multiple files are passed in, then they are
  processed in sync and if any name matches that read is discarded (or kept)
  from all files.
  
  positional arguments:
    fastqFiles            a (potentially gzipped) fastq file(s) containing the
                          reads with the order of reads the same in all files
  
  optional arguments:
    -h, --help            show this help message and exit
    -d DOTS, --dots DOTS  output dot to stderr every X reads. Input a negative
                          number to suppress output (default:-1)
    -o [OUTPUTFILES ...], --outputFiles [OUTPUTFILES ...]
                          an output file(s) (one for each input fastq file).
                          default(out1.fastq.gz ... outn.fastq.gz where n is the
                          number of fastqFiles)
    -f FILTERFILE, --filterFile FILTERFILE
                          a (potentially gzipped) file containing the names of
                          reads to be filtered one per line
    -k, --keep            keep reads matching the filter file and filter all
                          nonmatching reads
  
removeshort
~~~~

::
  
  usage: removeshort [-h] [-d DOTS] [-l MINLENGTH] [-n] [-p] [-b BADOUT]
                     fastqFile
  
  A program to remove short reads from a fastq file.
  
  positional arguments:
    fastqFile             a (potentially gzipped) fastq file containing the
                          sequence data
  
  optional arguments:
    -h, --help            show this help message and exit
    -d DOTS, --dots DOTS  output dot to stderr every X reads. Input a negative
                          number to suppress output (default:-1)
    -l MINLENGTH, --minLength MINLENGTH
                          minimum length read to output (default:15)
    -n, --removeN         remove reads which contain anything other than A, C, T
                          or G
    -p, --removePoor      remove reads with different length sequence and
                          qualities. Note this requires assuming that all reads
                          are 4 lines each
    -b BADOUT, --badOut BADOUT
                          a file path in which to save the first 10000 malformed
                          reads with different length sequence and qualities.
                          Note this requires assuming that all reads are 4 lines
                          each
  
splitbarcodes
~~~~

::
  
  usage: splitbarcodes [-h] [-i INDEXFILES [INDEXFILES ...]] [-d DOTS] -b
                       BARCODEFILE [-o OUTPUTPATH] [-u]
                       fastqFiles [fastqFiles ...]
  
  A program to take a list of barcodes and one or more fastq reads and one or
  two index reads and output reads matching the barcodes into a seperate file
  for each barcode. The script takes read files and index files where the reads
  and indexs are in the same order and outputs reads which match the appropriate
  barcodes into separate files.
  
  positional arguments:
    fastqFiles            a fastq file(s) (potentially gzipped) containing the
                          sequence reads
  
  optional arguments:
    -h, --help            show this help message and exit
    -i INDEXFILES [INDEXFILES ...], --indexFiles INDEXFILES [INDEXFILES ...]
                          a fastq file(s) (potentially gzipped) containing the
                          index reads
    -d DOTS, --dots DOTS  output dot to stderr every X reads. Input a negative
                          number to suppress output (default:-1)
    -b BARCODEFILE, --barcodeFile BARCODEFILE
                          a (potentially gzipped) file containing comma
                          separated sample names, first barcode and second
                          barcode (with no header and no commas in the sample
                          names)
    -o OUTPUTPATH, --outputPath OUTPUTPATH
                          a string giving the desired output directory
    -u, --unassigned      if set then store unassigned reads to
                          {outputPath}/__UNASSIGNED__R#.fastq.gz with their
                          corresponding barcodes in
                          {outputPath}/__UNASSIGNED__I#.fastq.gz
  
splitreadsbyname
~~~~

::
  
  usage: splitreadsbyname [-h] -s SPLITFILE [-d DOTS] [-o OUTPUTPATH] [-u] [-a]
                          fastqFile
  
  A program to take a list of read names and desired files and a fastq file
  containing various reads and output reads to their assigned file locations.
  The script takes a standard fastq read file and a csv separated file
  containing read identifiers and file locations and outputs reads which match
  the appropriate read names into their assigned files.
  
  positional arguments:
    fastqFile             a fastq file (potentially gzipped) containing the
                          sequence reads
  
  optional arguments:
    -h, --help            show this help message and exit
    -s SPLITFILE, --splitFile SPLITFILE
                          a (potentially gzipped) file containing comma
                          separated reads names and file names (with no header
                          and no commas in the sample names)
    -d DOTS, --dots DOTS  output dot to stderr every X reads. Input a negative
                          number to suppress output (default:-1)
    -o OUTPUTPATH, --outputPath OUTPUTPATH
                          a string giving the desired output directory
    -u, --unassigned      if set then store unassigned reads to
                          {outputPath}/__UNASSIGNED__.fastq.gz
    -a, --append          if set then append to already existing output files
  

Changelog
---------
0.1.7 (2023-02-11)

* Add `splitreadsbyname`

0.1.6 (2021-09-20)

* Fix code rot from python changes

0.1.5 (2018-12-21)

* Add `abitotrace` function

0.1.4 (2018-02-20)

* Adjust to changes in pysam v1.4
* Add option to only count bases with a quality greater above a specified limit

0.1.3 (2017-12-22)

* Add kmer counter
* Add option to filter reads containing N or unequal length sequence-qualities to removeshort

0.1.2 (2016-10-25)

* Add barcode splitter script

0.1.1 (2016-10-12)

* Add read filter script

0.1.0 (2016-01-20)

* Initial public release




