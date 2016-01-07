from setuptools import setup, find_packages
import sys
from setuptools.command.test import test as TestCommand

#https://pytest.org/latest/goodpractises.html
class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

setup(name='ampcountpy',
    version='0.1',
    description='Some handy functions for dealing with sequence alignment files.',
    url='http://github.com/sherrillmix/dnapy',
    author='Scott Sherrill-Mix',
    author_email='shescott@upenn.edu',
    license='GPL 3',
    packages=find_packages(),
    zip_safe=True,
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},
    entry_points={ 'console_scripts': [ 'ampcount = ampcountpy.__main__:main' ] },
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)'
    ],
    long_description=open('README.rst').read()
)
