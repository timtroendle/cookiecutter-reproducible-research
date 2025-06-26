![Reproduction](https://github.com/timtroendle/cookiecutter-reproducible-research/actions/workflows/reproduction.yaml/badge.svg)

# cookiecutter-reproducible-research

This repository provides [cookiecutter](http://cookiecutter.readthedocs.io) templates for reproducible research projects. The templates do not attempt to be generic, but have a clear and opinionated focus.

Projects build with these templates aim at full automation, and use `Python 3.11`, `conda`, `Git`, `Snakemake`, and `pandoc` to create a HTML and PDF report out of raw data, code, and `Markdown` text. Fork, clone, or download this repository on GitHub if you want to change any of these.

The template includes a few lines of code as a demo to allow you to create a report out of made-up simulation results right away. Read the `README.md` in the generated repository to see how.

These templates are developed on macOS and tested on Linux. They may work with Windows Subsystem for Linux, but Windows is not actively supported.

## Template types

> default

This generates the basic structure of a reproducible workflow.

> cluster

The cluster template extends the basic template by adding infrastructure to support running on a compute cluster.

## Getting Started

Make sure you have cookiecutter installed, otherwise install it with [conda](https://conda.io/docs/index.html):

    conda install cookiecutter -c conda-forge

Then create a repository using:

    cookiecutter gh:timtroendle/cookiecutter-reproducible-research --directory=[default/cluster]

You will be asked for the following parameters:

Parameter | Description
--- | ---
`project_name` | The name of your project, used in the documentation and report.
`project_short_name` | An abbreviation, used for environments and such. Avoid special characters and whitespace.
`author` | Your name.
`institute` | The name of your institute, used for report metadata.
`short_description` | A short description of the project, used for documentation and report.
`path_to_conda_envs` | The path to the directory hosting your conda envs (leave untouched for Snakemake default).

The `cluster` template requires the following parameter values in addition:

Parameter | Description
--- | ---
`cluster_url` | The address of the cluster to allow syncing to and from the cluster.
`cluster_base_dir` | The base path for the project on the cluster (default: `~/<project-short-name>`).
`cluster_type` | The type of job scheduler used on the cluster. Currently, only Slurm is supported.
`slurm_account` | The user account on Slurm.

## Project Structure

The generated repository will have the following structure:

```
├── config                  <- Configuration files, e.g., for your model if needed.
│   └── default.yaml        <- Default set of configuration parameter values.
├── data                    <- Raw input data.
├── envs                    <- Execution environments.
│   ├── default.yaml        <- Default execution environment.
│   ├── report.yaml         <- Environment for compilation of the report.
│   └── test.yaml           <- Environment for executing tests.
├── profiles                <- Snakemake profiles.
│   └── default             <- Default Snakemake profile folder.
│       └── config.yaml     <- Default Snakemake profile.
├── report                  <- All files creating the final report, usually text and figures.
│   ├── apa.csl             <- Citation style definition to be used in the report.
│   ├── literature.yaml     <- Bibliography file for the report.
│   ├── report.md           <- The report in Markdown.
│   └── pandoc-metadata.yaml<- Metadata for the report.
├── rules                   <- The place for all your Snakemake rules.
├── scripts                 <- Scripts go in here.
│   ├── model.py            <- Demo file.
│   └── vis.py              <- Demo file.
├── tests                   <- Automatic tests of the source code go in here.
│   └── test_model.py       <- Demo file.
├── .editorconfig           <- Editor agnostic configuration settings.
├── .ruff                   <- Linter and formatter settings for ruff.
├── .gitignore
├── environment.yaml        <- A file to create an environment to execute your project in.
├── LICENSE.md              <- MIT license description
├── Snakefile               <- Description of all computational steps to create results.
└── README.md
```

`cluster` templates additionally contain the following files:

```
├── envs
│   └── shell.yaml              <- An environment for shell rules.
├── profiles
│   └── cluster                 <- Cluster Snakemake profile folder.
│       └── config.yaml         <- Cluster Snakemake profile.
├── rules
│   └── sync.yaml               <- Snakemake rules to sync to and from the cluster.
├── .syncignore-receive         <- Build files to ignore when receiving from the cluster.
└── .syncignore-send            <- Local files to ignore when sending to the cluster.
```

## License

Some ideas for this cookiecutter template are taken from [cookiecutter-data-science](http://drivendata.github.io/cookiecutter-data-science/) and [mkrapp/cookiecutter-reproducible-science](https://github.com/mkrapp/cookiecutter-reproducible-science). This template is MIT licensed itself.
