"""This module contains the visualisation of results."""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def visualise_model_results(path_to_results, path_to_figure):
    """Plot the results."""
    sns.set_context('paper')
    results = pd.read_pickle(path_to_results)
    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    ax.plot(results)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.savefig(path_to_figure, dpi=300)


if __name__ == "__main__":
    visualise_model_results(
        path_to_results=snakemake.input.results[0],
        path_to_figure=snakemake.output[0]
    )
