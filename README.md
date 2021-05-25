![Reproduction](https://github.com/timtroendle/cookiecutter-reproducible-research/actions/workflows/reproduction.yaml/badge.svg)

# cookiecutter-reproducible-research

This repository provides a [cookiecutter](http://cookiecutter.readthedocs.io) template for reproducible research projects. It does not attempt to be generic, but has a clear and opinionated focus.

Projects build with this template aim at full automation, and use `Python 3.8`, `conda`, `Git`, `Snakemake`, and `pandoc` to create a HTML report out of raw data, code, and `Markdown` text. Fork, clone, or download this repository on GitHub if you want to change any of these.

The template includes a few lines of code as a demo to allow you to create a HTML report out of simulated results right away. Read the `README.md` in the generated repository to see how.

## Project Structure

The generated repository will have the following structure:

```
├── config                  <- Configuration files, e.g., for your model if needed.
├── data                    <- Raw input data.
├── envs                    <- Execution environments.
│   ├── default.yaml        <- Default execution environment.
│   ├── report.yaml         <- Environment for compilation of the report.
│   └── test.yaml           <- Environment for executing tests.
├── report                  <- All files creating the final report, usually text and figures.
│   ├── apa.csl             <- Citation style definition to be used in the report.
│   ├── literature.yaml     <- Bibliography file for the report.
│   ├── report.md           <- The report in Markdown.
│   └── pandoc-metadata.yaml<- Metadata for the report.
├── scripts                 <- Scripts go in here.
│   ├── model.py            <- Demo file.
│   └── vis.py              <- Demo file.
├── tests                   <- Automatic tests of the source code go in here.
│   └── test_model.py       <- Demo file.
├── .editorconfig           <- Editor agnostic configuration settings.
├── .flake8                 <- Linting settings for flake8.
├── .gitignore
├── environment.yaml        <- A file to create an environment to execute your project in.
├── LICENSE.md              <- MIT license description
├── Snakefile               <- Description of all computational steps to create results.
└── README.md
```

## Getting Started

Make sure you have cookiecutter installed, otherwise install it with [conda](https://conda.io/docs/index.html):

    conda install cookiecutter -c conda-forge

Then create a repository using:

    cookiecutter gh:timtroendle/cookiecutter-reproducible-research

You will be asked for the following parameters:

Parameter | Description
--- | ---
`project_name` | The name of your project, used in the documentation and report.
`project_short_name` | An abbreviation, used for environments and such. Avoid special characters and whitespace.
`author` | Your name.
`institute` | The name of your institute, used for report metadata.
`short_description` | A short description of the project, used for documentation and report.

## License

Some ideas for this cookiecutter template are taken from [cookiecutter-data-science](http://drivendata.github.io/cookiecutter-data-science/) and [mkrapp/cookiecutter-reproducible-science](https://github.com/mkrapp/cookiecutter-reproducible-science). This template is MIT licensed itself.
