# {{cookiecutter.project_name}}

{{cookiecutter.short_description}}

This repository contains the entire scientific project, including code and report. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

## Run the analysis

    snakemake --profile profiles/default

This will run all analysis steps to reproduce results and eventually build the report.

You can also run certain parts only by using other `snakemake` rules; to get a list of all rules run `snakemake --list`.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake --profile profiles/default -f dag

{% if cookiecutter._add_cluster_infrastructure == True -%}
## Run on a cluster

You may want to run the workflow on a cluster. While you can run on [any cluster that is supported by Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cluster.html), the workflow currently supports [LSF](https://en.wikipedia.org/wiki/Platform_LSF) clusters only. To run the workflow on a LSF cluster, use the following command:

    snakemake --profile profiles/cluster

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/cluster` as a starting point.

## Work local, build on remote

You may want to work locally (to change configuration parameters, add modules etc), but execute remotely on the cluster. This workflow supports you in working this way through three Snakemake rules: `send`, `receive`, and `clean_cluster_results`. It works like the following.

First, start local and make sure the `cluster-sync` configuration parameters fit your environment. Next, run `snakemake --profile profiles/default send` to send the entire repository to your cluster. On the cluster, execute the workflow with Snakemake (see above). After the workflow has finished, download results by locally running `snakemake --profile profiles/default receive`. By default, this will download results into `build/cluster`.

This workflow works iteratively too. After analysing your cluster results locally, you may want to make changes locally, send these changes to the cluster (`snakemake --profile profiles/default send`), rerun on the cluster, and download updated results (`snakemake --profile profiles/default receive`).

To remove cluster results on your local machine, run `snakemake --profile profiles/default clean_cluster_results`.
{%- endif %}

## Be notified of build successes or fails

  As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully. Notifications are sent by email. To activate notifications, add the email address of the recipient to the configuration key `email`. You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

      snakemake --profile profiles/default --config email=<your-email>

## Run the tests

    snakemake --profile profiles/default test

## Repo structure

* `report`: contains all files necessary to build the report; plots and result files are not in here but generated automatically
* `scripts`: contains the Python source code as scripts
* `rules`: contains Snakemake rule definitions
* `envs`: contains execution environments
* `tests`: contains the test code
* `config`: configurations used in the study
* `profiles`: Snakemake execution profiles
* `data`: place for raw data
* `build`: will contain all results (does not exist initially)

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
