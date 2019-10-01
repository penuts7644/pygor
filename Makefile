##
##	ImmunoProbs
##		Create IGoR models and calculate the generation probability of V(D)J
##		and CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt
##

VERSION = $(shell git tag | tail -1)

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
	pip install -I -r requirements_travis.txt && pip install -e .[development]

##		make test
##			Run pytest tests from the tests directory on the immuno_probs source.
##
test:
	python -m pytest -v tests

##		make clean
##			Removes the old distribution directories and files.
##
clean:
	rm -rf ./dist && rm -rf ./build && rm -rf ./docs/_build && rm -rf ./.pytest_cache && find . -name '*.pyc' -type f -delete

##		make build
##			Perfoms tests, a dir clean, builds the new distribution package as
##			well as the documentation.
##
build: test clean
	python setup.py bdist_wheel

##		make docs
##			Build the documentation for ImmunoProbs.
##
docs: test clean
	cd docs && make html

##		make build-docker
##			Perfoms tests, a dir clean, builds the new distribution package for
##			ImmunoProbs and finally builds a docker image of all executables.
##
build-docker: test clean build
	docker build -t penuts7644/immuno-probs:$(VERSION) . && docker tag penuts7644/immuno-probs:$(VERSION) penuts7644/immuno-probs:latest

##		make test-deploy
##			Tests, cleans, builds and uploads all distribution files to PyPI test server.
##
test-deploy: test clean build
	python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*

##		make deploy
##			Tests, cleans, builds and uploads all distribution files to PyPI.
##
deploy: test clean build
	python -m twine upload dist/*
