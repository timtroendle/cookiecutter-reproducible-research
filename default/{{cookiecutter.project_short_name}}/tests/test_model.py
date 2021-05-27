"""Test case for the model."""
import scripts.model


def test_model():
    assert scripts.model.linear_model(slope=1, x0=4, x=0) == 4
