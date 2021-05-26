"""This package contains the model used in this study."""
import numpy as np
import pandas as pd


def linear_model(slope, x0, x):
    """Returns the function value of a linear model.

    The model is y=slope*x + x0.

    Parameters:
        * slope: the slope of the linear function
        * x0: the y-value at x=0
        * x: the x value for which the y value should be determined
    """
    return slope * x + x0


def computation_step(slope, x0, path_to_results):
    """Evaluates the linear model between x=0 and x=1 and writes results."""
    x = np.linspace(start=0, stop=1, num=50)
    result = linear_model(slope, x0, x)
    pd.Series(index=x, data=result).to_pickle(path_to_results)


if __name__ == "__main__":
    computation_step(
        slope=snakemake.params.slope,
        x0=snakemake.params.x0,
        path_to_results=snakemake.output[0]
    )
