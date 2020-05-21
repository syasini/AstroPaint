"""This module contains shared test fixtures such as Catalog and Canvas instances"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import pytest
from astropaint import Catalog, Canvas


@pytest.fixture(scope="session")
def nside():
    yield 64


@pytest.fixture(scope="module")
def test_catalog():
    """Create a test catalog to be used in unit tests"""
    print(f"\n{'':->10}initializing the catalog{'':-<10}\n")

    catalog = Catalog("test")
    yield catalog

    print(f"\n{'':->10} closing the catalog{'':-<10}\n")


@pytest.fixture(scope="function")
def test_canvas(test_catalog, nside):
    """Create a test canvas to be used in unit tests"""
    print(f"\n{'':->10}initializing the canvas{'':-<10}\n")

    canvas = Canvas(test_catalog, nside)
    yield canvas

    print(f"\n{'':->10} closing the canvas{'':-<10}\n")