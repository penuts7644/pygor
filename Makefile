##
##	pygor
##		Pygor is part of the IGoR (Inference and Generation of Repertoires) software.
##		Pygor Python package can be used to post process files generated by IGoR.
##		Copyright (C) 2018 Quentin Marcou & Wout van Helvoirt
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
##			requirements.
##
setup:
	pip install -e .[development]

##		make pytest
##			Run pytest tests from the tests directory on the pygor source.
##
pytest:
	python -m pytest -v tests

##		make test-build
##			Perfoms clean and builds the new distribution package.
##
test-build:
	rm -rf ./dist && rm -rf ./build && find . -name '*.pyc' -type f -delete
	python setup.py sdist bdist_wheel

##		make test-deploy
##			Upload all distribution files to PyPI test server.
##
test-deploy:
	python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
