"""This module contains unit tests for the Catalog class"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import pytest


def test_alm_Cl_outdated_upon_clean_canvas(test_canvas):
    """Passes if after clean.canvas()
    the status of pixels, alms, and Cls are correctly updated """
    test_canvas.clean()
    assert not test_canvas._pixel_is_outdated
    assert test_canvas._alm_is_outdated
    assert test_canvas._Cl_is_outdated


def test_alm_Cl_outdated_upon_pixel_change(test_canvas):
    """Passes if after setting canvas.pixels
    the status of pixels, alms, and Cls are correctly updated """
    test_canvas.pixels = 1
    assert not test_canvas._pixel_is_outdated
    assert test_canvas._alm_is_outdated
    assert test_canvas._Cl_is_outdated


def test_pixel_Cl_outdated_upon_alm_change(test_canvas):
    """Passes if after setting canvas.alms
    the status of pixels, alms, and Cls are correctly updated """
    test_canvas.alm = 1
    assert test_canvas._pixel_is_outdated
    assert not test_canvas._alm_is_outdated
    assert test_canvas._Cl_is_outdated


def test_error_upon_Cl_change(test_canvas):
    """Passes if upon setting Cl TypeError is raised"""
    try:
        test_canvas.Cl = 1
    except AttributeError:  # Cl is a property
        assert True

@pytest.mark.parametrize("resolution", [2, 1/2])
def test_alm_Cl_outdated_upon_ud_grade(test_canvas, nside, resolution):
    """Passes if after upgrading and downgrading resolution
    the status of pixels, alms, and Cls and lmax are correctly updated """
    print(test_canvas.alm)
    # test various resolution
    new_nside = int(nside*resolution)
    test_canvas.ud_grade(new_nside)

    assert not test_canvas._pixel_is_outdated
    assert test_canvas._alm_is_outdated
    assert test_canvas._Cl_is_outdated
    assert test_canvas.lmax == 3 * new_nside - 1

