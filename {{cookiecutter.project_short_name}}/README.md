# {{cookiecutter.project_name}}

{{cookiecutter.short_description}}

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

The following dependencies are needed to set up an environment in which the analysis can be run and the paper be build:

* [conda](https://conda.io/docs/index.html)
* `LaTeX` to [produce a PDF](http://pandoc.org/MANUAL.html#creating-a-pdf). Can be avoided by switching to [any other output format supported by pandoc](http://pandoc.org/index.html).
* `make` (optional; without `make`, e.g. on Windows, you will need to manually run all steps)

When these dependencies are installed, you can create a conda environment from within you can run the analysis:

    conda env create -f conda-environment.yml

Don't forget to activate the environment. To see what you can do now, run:

    make help

## Run the analysis

    make paper

This will run all analysis steps to reproduce results and eventually build the paper.

You can also run certain parts only by using other `make` rules; to get a list of all rules see the `Makefile`.

If you do not have `make` you can manually run the steps through the Python command line interfaces. Refer to the `Makefile` to see which commands are called to produce results.

## Run the tests

    make test

## Repo structure

* `report`: contains all files necessary to build the paper; plots and result files are not in here but generated automatically
* `src`: contains the Python source code
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)
