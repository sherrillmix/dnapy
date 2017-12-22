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

INSERT_USAGE_HERE

Changelog
---------
0.1.3 (2017-12-22)

* Add kmer counter
* Add option to filter reads containing N or unequal length sequence-qualities to removeshort

0.1.2 (2016-10-25)

* Add barcode splitter script

0.1.1 (2016-10-12)

* Add read filter script

0.1.0 (2016-01-20)

* Initial public release




