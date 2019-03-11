##
##	ImmunoProbs
##		ImmunoProbs Python package able to calculate the generation probability
##		of V(D)J and CDR3 sequences. Copyright (C) 2019 Wout van Helvoirt
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

##		make build-docker
##			Perfoms tests, a dir clean, builds the new distribution package for
##			ImmunoProbs and finally builds a docker image of all executables.
##
build-docker: test clean build
	docker build -t penuts7644/immuno-probs:0.1.6 . && docker tag penuts7644/immuno-probs:0.1.6 penuts7644/immuno-probs:latest

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
