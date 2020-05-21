"""This module contains unit tests for the Catalog class"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"


def test_correct_input(test_catalog):
    """Passes if coordinates, velocities and Mass is in the input catalog"""
    attributes = ["x", "y", "z", "v_x", "v_y", "v_z", "M_200c"]

    assert all([item in test_catalog.data for item in attributes ])


def test_catalog_size(test_catalog):
    """Passes if catalog.size matches the length of data"""

    assert test_catalog.size == len(test_catalog.data)


def test_has_data(test_catalog):
    """Passes if catalog has a data attribute"""
    try:
        test_catalog.data
    except AttributeError:
        assert False

