##
##	ImmunoProbs
##		ImmunoProbs Python package uses a simplified manner for calculating the
##		generation probability of V(D)J and CDR3 sequences.
##		Copyright (C) 2018 Wout van Helvoirt
##

default: help

##	COMMANDS
##

##		make help
##			Display the help information.
##
help:
	@grep '^##.*' ./Makefile

##		make setup
##			Setup the development enviroment from setup.py and install all
##			additional development requirements.
##
setup:
	pip install -r requirements_travis.txt && pip install -e .[development]

##		make test
##			Run pytest tests from the tests directory on the immuno_probs source.
##
test:
	python -m pytest -v tests

##		make clean
##			Removes the old distribution directories and files.
##
clean:
	rm -rf ./dist && rm -rf ./build && find . -name '*.pyc' -type f -delete

##		make build
##			Perfoms tests, a dir clean and builds the new distribution package.
##
build: test clean
	python setup.py bdist_wheel

##		make test-deploy
##			Upload all distribution files to PyPI test server.
##
test-deploy:
	python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

##		make deploy
##			Tests, cleans, builds and upload all distribution files to PyPI.
##
deploy: test clean build
	python -m twine upload dist/*
