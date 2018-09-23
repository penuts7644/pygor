# Pygor is part of the IGoR (Inference and Generation of Repertoires)
# software. This Python package can be used for post processing of IGoR
# generated output files.
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


from setuptools import setup

version_num = "0.1.0"

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("LICENSE") as f:
    license = f.read()

with open("requirements.txt") as f:
    requirements = [line.rstrip("\n") for line in f.readlines()]

setup(
    name="pygor",
    version=version_num,
    description="This Python package can be used for post "
                "processing of IGoR generated output files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/penuts7644/pygor",
    download_url="https://github.com/penuts7644/pygor/archive/v"
                 + version_num + ".tar.gz",
    author="Wout van Helvoirt",
    author_email="wout.van.helvoirt@icloud.com",
    license=license,
    keywords=["IGoR", "Bio-Informatics", "V(D)J Recombination", "Sequencing",
              "Analysis", "DNA", "Models"],
    packages=[
        "pygor",
        "pygor.aligns",
        "pygor.counters",
        "pygor.counters.bestscenarios",
        "pygor.counters.coverage",
        "pygor.models",
        "pygor.utils"
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],
    install_requires=requirements,
    zip_safe=False)
