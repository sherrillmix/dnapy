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

Usage
-----
The package provides the scripts:

splitbarcodes
~~~~

::
  
  usage: splitbarcodes [-h] [-i INDEXFILES [INDEXFILES ...]] [-d DOTS] -b
                       BARCODEFILE [-o OUTPUTPATH]
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
                          a file (potentially gzipped) file containing comma
                          separated sample names, first barcode and second
                          barcode (with no header and no commas in the sample
                          names)
    -o OUTPUTPATH, --outputPath OUTPUTPATH
                          a string giving the desired output directory
  
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
  
countbases
~~~~

::
  
  usage: countbases [-h] [-v] [-r REGION] [-s] bamFile
  
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
  
removereads
~~~~

::
  
  usage: removereads [-h] [-d DOTS] -f FILTERFILE [-o OUTPUTFILES]
                     fastqFiles [fastqFiles ...]
  
  A program to filter reads by name from a single/set of fastq file(s). The
  script looks for reads which have a name line where the string before a space
  exactly matches a pattern. If multiple files are passed in, then they are
  processed in sync and if any name matches that read is discarded from all
  files.
  
  positional arguments:
    fastqFiles            a fastq or fastqs (potentially gzipped) file
                          containing the reads
  
  optional arguments:
    -h, --help            show this help message and exit
    -d DOTS, --dots DOTS  output dot to stderr every X reads. Input a negative
                          number to suppress output (default:-1)
    -f FILTERFILE, --filterFile FILTERFILE
                          a file (potentially gzipped) file containing the names
                          of reads to be filtered one per line
    -o OUTPUTFILES, --outputFiles OUTPUTFILES
                          a comma separated string giving a list of output files
                          (one for each input fastq file). default(out1.fastq.gz
                          ... outn.fastq.gz where n is the number of fastqFiles)
  
removeshort
~~~~

::
  
  usage: removeshort [-h] [-d DOTS] [-l MINLENGTH] fastqFile
  
  A program to remove short reads from a fastq file.
  
  positional arguments:
    fastqFile             a fastq (potentially gzipped) file containing the
                          alignment
  
  optional arguments:
    -h, --help            show this help message and exit
    -d DOTS, --dots DOTS  output dot to stderr every X reads. Input a negative
                          number to suppress output (default:-1)
    -l MINLENGTH, --minLength MINLENGTH
                          minimum length read to output (default:15)
  

Changelog
---------
0.1.1 (2016-10-12)

* Add read filter script

0.1.0 (2016-01-20)

* Initial public release




