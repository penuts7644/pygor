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

##		make test
##			Test immuno_probs source code with flake8 linter and pytest.
##
test:
	flake8 . --count --select=E9,F63,F7,F82 --show-source --statistics && flake8 . --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics && python -m pytest -v tests

##		make tag v=<*.*.*>
##			Generate changelog file, tag the latest commit with specified version
##			 and push the tag.
##
tag:
	if [ -z "$v" ]; then echo "Argument missing or empty: 'v=*.*.*'"; else git log $(VERSION)..HEAD --pretty=format:"%s" -i -E > change-log.txt && git tag $(v) && git push --tags; fi

##		make build-docs
##			Build the documentation for ImmunoProbs.
##
build-docs: test clean
	cd docs && make html

##		make build-pypi
##			Perfom tests, a directory cleanup and a build of the ImmunoProbs python
##			distribution package.
##
build-pypi: test clean
	python setup.py bdist_wheel

##		make build-docker
##			Perfom tests, a directory cleanup, build the python distribution package
##			for ImmunoProbs and build a docker image of all executables.
##
build-docker: build-pypi
	docker build -t penuts7644/immuno-probs:$(VERSION) . && docker tag penuts7644/immuno-probs:$(VERSION) penuts7644/immuno-probs:latest

##		make deploy-pypi
##			Deploy ImmunoProbs python distribution files to PyPI.
##
deploy-pypi:
	python -m twine upload dist/*

##		make deploy-docker
##			Deploy the ImmunoProbs docker images (<version> and 'latest') to
##			docker hub.
##
deploy-docker:
	docker push penuts7644/immuno-probs:$(VERSION) && docker push penuts7644/immuno-probs:latest
