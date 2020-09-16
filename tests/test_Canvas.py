"""This module contains unit tests for the Catalog class"""

__author__ = "Siavash Yasini"
__email__ = "yasini@usc.edu"

import pytest
from .conftest import apply_func_to_cutout

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

def test_pixels_are_zero_after_clean(test_canvas):
    """Passes if all pixels are zero after canvas.clean()"""
    # suggested by dsspiegel:
    test_canvas.clean()
    assert (test_canvas.pixels == 0).all(), f"Pixels not all zero: {test_canvas.pixels}"


@pytest.mark.parametrize("mult_by, add_to",
                         [(1, 1),  # both scalars
                          (1, [1]),  # one scalar one vector
                          ([1], [1]),  # two vectors
                          (1, [1, 2, 3, 4, 5, 6]),  # vectors of different sizes
                          ])
def test_cutouts_apply_func_kwarg_types(test_canvas, mult_by, add_to):
    """Passes if various kwarg types (scalar, list, etc) can be passed to apply func"""

    # get the halo list from catalog
    halo_list = test_canvas.catalog.data.index

    # cutout patches and apply test function to them
    cutouts = test_canvas.cutouts(halo_list,
                                  apply_func=apply_func_to_cutout,
                                  func_kwargs={"mult_by": mult_by,
                                               "add_to": add_to})
    # iterate through the cutouts
    for _ in cutouts:
        pass


def test_cutouts_apply_multiple_funcs(test_canvas):
    """Passes if multiple apply funcs can be applied to patch"""

    # get the halo list from catalog
    halo_list = test_canvas.catalog.data.index

    # cutout patches and apply test function to them
    cutouts = test_canvas.cutouts(halo_list,
                                  apply_func=[apply_func_to_cutout]*2,
                                  func_kwargs=[{"mult_by": 1,
                                               "add_to": 1}]*2)
    # iterate through the cutouts
    for _ in cutouts:
        pass
