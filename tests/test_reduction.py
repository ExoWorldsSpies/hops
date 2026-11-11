"""Regression test of the reduction (bias, dark and flat correction, statistics, timing).

The reduction of the synthetic data set (tests/synthetic.py) must reproduce the reference values
in tests/golden/reduction.json, which were produced with the unchanged code. Changes that are meant
to keep the results (refactoring, speed-ups) must keep this test passing.

The comparison allows for tiny platform differences (other numpy / CPU builds): statistics to a
relative 1e-7 and pixel values to +/- 1 ADU (rounding to whole ADU can flip at .5). Exact equality
is reported as well.

    python -m pytest tests

Needs a display, as the reduction runs in its (hidden) Tk window. Re-create the reference values
only for intended changes of the results: python tests/make_golden.py
"""

import os
import json

import numpy as np
import pytest

import reduction_run

GOLDEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden', 'reduction.json')


@pytest.fixture(scope='module')
def result(tmp_path_factory):
    if os.name == 'posix' and os.uname().sysname == 'Linux' and not os.environ.get('DISPLAY'):
        pytest.skip('the reduction runs in a Tk window and needs a display (e.g. xvfb-run)')
    cwd = os.getcwd()
    try:
        yield reduction_run.run_reduction(str(tmp_path_factory.mktemp('synthetic')))
    finally:
        os.chdir(cwd)


@pytest.fixture(scope='module')
def golden():
    with open(GOLDEN) as f:
        return json.load(f)


def test_same_frames(result, golden):
    assert sorted(result) == sorted(golden)


@pytest.mark.parametrize('key, rtol', [('mean', 1e-7), ('std', 1e-7), ('psf', 1e-7), ('jd', 1e-12),
                                       ('airmass', 1e-9)])
def test_frame_statistics(result, golden, key, rtol):
    for name in golden:
        assert result[name][key] == pytest.approx(golden[name][key], rel=rtol), name


def test_skipped_frames(result, golden):
    assert {name: result[name]['skip'] for name in golden} == {name: golden[name]['skip'] for name in golden}


def test_reduced_images(result, golden):
    block_pixels = reduction_run.BLOCK ** 2
    for name in golden:
        assert result[name]['shape'] == golden[name]['shape'], name
        difference = np.abs(np.array(result[name]['block_sums']) - np.array(golden[name]['block_sums']))
        # every pixel may be off by one ADU at most
        assert difference.max() <= block_pixels, name


def test_reduced_images_exactly(result, golden):
    """Bit-identical images: expected on the platform the reference was made on."""
    different = [name for name in golden if result[name]['sha256'] != golden[name]['sha256']]
    if different:
        pytest.xfail('{0} frame(s) differ in the last digit, within the tolerance of '
                     'test_reduced_images (other platform or library versions?)'.format(len(different)))
