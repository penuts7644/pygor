##
##	ImmunoProbs
##		Create IGoR models and calculate the generation probability of V(D)J and CDR3 sequences.
##		Copyright (C) 2019 Wout van Helvoirt
##

VERSION = $(shell git tag | tail -1)

default: help

##	COMMANDS
##

##		make help
##			Display the help information for this file.
##
help:
	@grep '^##.*' ./Makefile

##		make clean
##			Delete the old distribution and build artifacts.
##
clean:
	rm -rf ./dist && rm -rf ./build && rm -rf ./docs/_build && rm -rf ./.pytest_cache && find . -name '*.pyc' -type f -delete

##		make setup
##			Perform directory cleanup and setup the development enviroment.
##
setup: clean
	pip install -I -r requirements.txt && pip install -e .[development]

##		make tag v=<v*.*.*>
##			Generate changelog file, tag the latest commit with specified version and push the tag.
##
tag:
	if [ -z "$v" ]; then echo "Argument missing or empty: 'v=v*.*.*'"; else git log $(VERSION)..HEAD --pretty=format:"%s" -i -E > change-log.txt && git tag $(v) && git push --tags; fi

##		make docs
##			Build the documentation for ImmunoProbs.
##
docs: clean
	cd docs && make html

##		make build
##			Perfom tests, a directory cleanup and a build of the ImmunoProbs python	distribution package.
##
build: clean
	python setup.py bdist_wheel

##		make deploy
##			Deploy ImmunoProbs python distribution files to PyPI.
##
deploy:
	python -m twine upload dist/*
