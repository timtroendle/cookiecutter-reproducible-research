# {{cookiecutter.project_name}}

{{cookiecutter.short_description}}

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

The following dependencies are needed to set up an environment in which the analysis can be run and the report be build:

* [conda](https://conda.io/docs/index.html)
* `LaTeX` to [produce a PDF](http://pandoc.org/MANUAL.html#creating-a-pdf). Can be avoided by switching to [any other output format supported by pandoc](http://pandoc.org/index.html).

When these dependencies are installed, you can create a conda environment from within you can run the analysis:

    conda env create -f conda-environment.yml

Don't forget to activate the environment. To see what you can do now, run:

    snakemake --list

## Run the analysis

    snakemake report

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps, and if you have `dot` installed, run:

    snakemake --rulegraph | dot -Tpdf > dag.pdf

## Run the tests

    snakemake test

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `src`: contains the Python source code
* `tests`: contains the test code
* `config`: configurations used in the study
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)
