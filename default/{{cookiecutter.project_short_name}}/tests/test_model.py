"""Test case for the model results.

Pytest fixtures are provided from within the test runner.
"""


def test_model_results_at_0(model_results):
    assert model_results[0] == 5


def test_model_results_at_1(model_results):
    assert model_results[1] == 9
