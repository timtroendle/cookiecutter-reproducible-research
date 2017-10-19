# cookiecutter-reproducible-science

This repository provides a [cookiecutter](http://cookiecutter.readthedocs.io) template for reproducible science projects. It does not attempt to be generic, but has a clear and opinionated focus.

Projects build with this template aim at full automation, and use `Python 3.6`, `conda`, `Git`, `Make`, `pandoc`, and `LaTeX` to create a PDF report out of raw data, code, and `Markdown` text. Fork, clone, or download this repository on GitHub if you want to change any of these.

The template includes a few lines of code as a demo to allow you to create a PDF report out of simulated results right away. Read the `README.md` in the generated repository to see how.

## Project Structure

The generated repository will have the following structure:

```
├── config                  <- Configuration files, e.g., for your model if needed.
├── data                    <- Raw input data.
├── report                  <- All files creating the final report, usually text and figures.
│   ├── energy-policy.csl   <- Citation style definition to be used in the report.
│   ├── literature.bib      <- Bibliography file for the report.
│   ├── main.md             <- The report in Markdown.
│   └── pandoc-metadata.yml <- Metadata for the report.
├── src                     <- Source code goes in here.
│   ├── __init__.py         <- Makes `src` a Python module.
│   ├── model.py            <- Demo file.
│   └── vis.py              <- Demo file.
├── tests                   <- Automatic tests of the source code go in here.
│   └── test_model.py       <- Demo file.
├── .gitignore
├── conda-environment.yml   <- A file to create an environment to execute your project in.
├── LICENSE.md              <- MIT license description
├── Makefile                <- Description of all computational steps to create results.
├── README.md
└── VERSION
```

## Getting Started

Make sure you have cookiecutter installed, otherwise install it with [conda](https://conda.io/docs/index.html):

    conda install cookiecutter -c conda-forge

Then create a repository using:

    cookiecutter gh:timtroendle/cookiecutter-reproducible-science

You will be asked for the following parameters:

Parameter | Description
--- | ---
`project_name` | The name of your project, used in the documentation and report.
`project_short_name` | An abbreviation, used for environments and such. Avoid special characters and whitespace.
`author` | Your name.
`institute` | The name of your institute, used for report metadata.
`short_description` | A short description of the project, used for documentation and report.
`version` | The version of the project.

## License

Some ideas for this cookiecutter template are taken from [cookiecutter-data-science](http://drivendata.github.io/cookiecutter-data-science/) and [mkrapp/cookiecutter-reproducible-science](https://github.com/mkrapp/cookiecutter-reproducible-science). This template is MIT licensed itself.
