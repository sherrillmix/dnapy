ampCountPy
==========
.. image:: https://travis-ci.org/sherrillmix/ampCountPy.svg?branch=master
    :alt: Travis CI
    :target: https://travis-ci.org/sherrillmix/ampCountPy
.. image:: https://codecov.io/github/sherrillmix/ampCountPy/coverage.svg?branch=master
    :alt: Codecov
    :target: https://codecov.io/github/sherrillmix/ampCountPy?branch=master
.. image:: https://badge.fury.io/py/ampcountpy.svg
    :alt: PyPi version
    :target: https://badge.fury.io/py/ampcountpy


Some python functions to count the expected amplifications for genomic regions given a set of primer binding locations for a `multiple displacement amplification <http://en.wikipedia.org/wiki/Multiple_displacement_amplification>`_ reaction. See `ampCountR <https://github.com/sherrillmix/ampCountR>`_ for more details.
 
Installation
------------
Easy install
~~~~~~~~~~~~

The easy way to install is to just do::

  pip install ampcountpy

Github
~~~~~~

To install the development version from github, clone the repository to a local directory using something like::

    git clone https://github.com/sherrillmix/ampcountpy.git

and run `setup.py` from the resulting directory (the `--user` installs it locally and doesn't require root access)::

  cd ampcountpy
  python setup.py install --user
  python setup.py test

Run directly
------------
The module can be called directly using something like::

  python -m ampcountpy -f forward.txt -r reverse.txt

or::

  ampcountpy -f forward.txt -r reverse.txt

where `forward.txt` is a text file containing position of primer landing sites on the forward strand and `reverse.txt` is primer landing sites on the reverse strand. By default, amplification predictions are output to out.csv. The full details on options and arguments is available with::

  ampcountpy --help

Using function in python
------------------------
The main function is `predictAmplifications` which can be used like:

..  code:: python

  from ampcountpy import predictAmplifications
    forwards=[1,2,3]
    reverses=[5,6,7]
    predictions=predictAmplifications(forwards,reverses)

where `forwards` are the 5'-most base of primer landing sites on the forward strand and `reverses` are the 3'-most base of primers landing on the reverse strand.


Changelog
---------
0.1.3 (2015-11-02)

* Fix header

0.1.2 (2015-11-02)

* Fix changelog formatting

0.1.1 (2015-11-02)

* Pip install instructions

0.1.0 (2015-11-02)

* Initial public release




