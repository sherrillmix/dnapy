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

INSERT_USAGE_HERE

Changelog
---------
0.1.9 (2022-03-20)

* Add helper class `Trie`

0.1.8 (2023-03-06)

* Add `countbarcodes`

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




