"""Test case for the model."""
import src.model


def test_model():
    assert src.model.linear_model(slope=1, x0=4, x=0) == 4
