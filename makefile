.PHONY: install test

README.rst: dnapy/*.py preREADME.rst setup.py makeReadme.py
	#need a readme to install
	cp preREADME.rst README.rst
	make install
	python makeReadme.py
	make install

test:
	python setup.py test
install:
	python setup.py install --user
