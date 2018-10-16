# Pygor is part of the IGoR (Inference and Generation of Repertoires) software.
# Pygor Python package can be used to post process files generated by IGoR.
# Copyright (C) 2018 Quentin Marcou & Wout van Helvoirt

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from setuptools import setup, find_packages


VERSION = '1.0.0-test'

with open('README.md', 'r') as f:
    LONG_DESCRIPTION = f.read()

with open('requirements.txt') as f:
    REQUIREMENTS = [line.rstrip('\n') for line in f.readlines()]

setup(
    name='pygor',
    version=VERSION,
    description='Pygor Python package can be used to post process files generated by IGoR.',
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url='https://github.com/penuts7644/pygor',
    author='Wout van Helvoirt',
    author_email='wout.van.helvoirt@icloud.com',
    keywords='IGoR Bio-Informatics Recombination Sequencing Analysis DNA Models Genes',
    packages=find_packages(include=["pygor", "pygor.*"]),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
    ],
    install_requires=REQUIREMENTS,
    zip_safe=False,
    project_urls={
        'IGoR Source': 'https://github.com/qmarcou/IGoR',
        'Bug Reports': 'https://github.com/penuts7644/pygor/issues',
        'Release Notes': 'https://github.com/penuts7644/pygor/releases/tag/v' + VERSION,
    },
)
