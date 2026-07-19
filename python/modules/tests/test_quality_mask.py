"""Regression tests for :func:`emit_tools.quality_mask`.

These tests are network-free: they build small synthetic EMIT L2A Mask netCDF
files whose ``sensor_band_parameters/mask_bands`` metadata mirrors the real
V001 (8-band) and V002 (11-band) products. The V002 band names and ordering
were verified against a real granule
(``EMIT_L2A_MASK_002_20220810T034103_2222203_001.nc``):

    0  Cloud Flag                 (binary flag)
    1  Cirrus Flag                (binary flag)
    2  Water Flag                 (binary flag)
    3  Spacecraft Flag            (binary flag)
    4  Dilated Cloud Flag         (binary flag)
    5  AOD550                     (continuous data)
    6  H2O (g cm-2)               (continuous data)
    7  Aggregate Flag             (binary flag)
    8  SpecTf-Cloud Probability   (continuous data, probability)
    9  SpecTf-Cloud Flag          (binary flag)
    10 SpecTf-Buffer Distance     (continuous data, distance)

The bug being guarded against: the previous implementation hard-coded bands 5
and 6 as the only non-flag layers, so it silently accepted the new V002
continuous bands (8 and 10) and returned a semantically meaningless "mask".

An optional test also runs against a real granule when the environment variable
``EMIT_L2A_MASK_V002`` points to one.
"""

import os

import numpy as np
import netCDF4 as nc
import pytest

from emit_tools import quality_mask


# --- Verified real band layouts -------------------------------------------------

V002_MASK_BANDS = [
    "Cloud Flag",
    "Cirrus Flag",
    "Water Flag",
    "Spacecraft Flag",
    "Dilated Cloud Flag",
    "AOD550",
    "H2O (g cm-2)",
    "Aggregate Flag",
    "SpecTf-Cloud Probability",
    "SpecTf-Cloud Flag",
    "SpecTf-Buffer Distance",
]

# Representative EMIT L2A Mask V001 layout (8 bands). V001 mixed the
# capitalization of "flag"/"Flag", which the classifier must handle
# case-insensitively. AOD550/H2O sit at indices 5/6, matching the constants the
# previous implementation hard-coded.
V001_MASK_BANDS = [
    "Cloud flag",
    "Dilated Cloud Flag",
    "Cirrus flag",
    "Water flag",
    "Spacecraft Flag",
    "AOD550",
    "H2O (g cm-2)",
    "Aggregate Flag",
]


def _write_mask_nc(path, band_names, mask_array):
    """Write a minimal EMIT-like L2A Mask netCDF file.

    Structure mirrors the real product: a root ``mask`` variable with dims
    ``(downtrack, crosstrack, bands)`` and a ``sensor_band_parameters`` group
    holding the ``mask_bands`` string variable.
    """
    downtrack, crosstrack, n_bands = mask_array.shape
    assert n_bands == len(band_names)
    with nc.Dataset(path, "w", format="NETCDF4") as ds:
        ds.createDimension("downtrack", downtrack)
        ds.createDimension("crosstrack", crosstrack)
        ds.createDimension("bands", n_bands)
        mask_var = ds.createVariable(
            "mask", "f4", ("downtrack", "crosstrack", "bands")
        )
        mask_var[:] = mask_array.astype("f4")
        mask_var.long_name = "Masks"
        mask_var.units = "unitless"
        grp = ds.createGroup("sensor_band_parameters")
        band_var = grp.createVariable("mask_bands", str, ("bands",))
        for i, name in enumerate(band_names):
            band_var[i] = name
        band_var.long_name = "Mask Band Names"
    return str(path)


def _binary(pattern):
    return np.array(pattern, dtype="f4")


# Reusable binary flag patterns on a 3x4 grid.
_A = _binary([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]])
_B = _binary([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
_C = _binary([[1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]])
_D = _binary([[1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
_E = _binary([[0, 1, 0, 1], [0, 1, 1, 0], [1, 0, 0, 0]])
_ZERO = np.zeros((3, 4), dtype="f4")

# Continuous layers (deliberately include values > 1 and fractional values).
_AOD = _binary([[-0.04, 0.20, 0.50, 1.10], [0.30, 0.80, 0.95, 1.05], [0.00, 0.10, 0.90, 0.99]])
_H2O = _binary([[0.05, 1.20, 2.00, 3.60], [0.50, 1.10, 0.90, 2.20], [0.30, 0.40, 1.50, 0.80]])
_PROB = _binary([[0.30, 0.95, 0.50, 1.00], [0.26, 0.70, 0.99, 0.40], [0.80, 0.35, 0.60, 0.45]])
# Buffer distance is conceptually continuous (pixel distance); include 2 and 3.
_DIST = _binary([[0, 1, 2, 3], [0, 0, 1, 2], [3, 2, 1, 0]])


@pytest.fixture
def v002_mask(tmp_path):
    layers = [_A, _B, _ZERO, _ZERO, _C, _AOD, _H2O, _D, _PROB, _E, _DIST]
    arr = np.stack(layers, axis=-1)
    path = _write_mask_nc(tmp_path / "v002_mask.nc", V002_MASK_BANDS, arr)
    return {"path": path, "names": V002_MASK_BANDS, "array": arr}


@pytest.fixture
def v001_mask(tmp_path):
    layers = [_A, _B, _C, _ZERO, _ZERO, _AOD, _H2O, _D]
    arr = np.stack(layers, axis=-1)
    path = _write_mask_nc(tmp_path / "v001_mask.nc", V001_MASK_BANDS, arr)
    return {"path": path, "names": V001_MASK_BANDS, "array": arr}


def _expected_or(arr, bands):
    """Reference mask: logical OR of the given (binary) band layers."""
    out = np.zeros(arr.shape[:2], dtype=np.uint8)
    for b in bands:
        out |= (arr[:, :, b] > 0).astype(np.uint8)
    return out


def _assert_binary_mask(qmask):
    assert isinstance(qmask, np.ndarray)
    assert qmask.dtype == np.uint8
    assert set(int(v) for v in np.unique(qmask)).issubset({0, 1})


# --- V001: legacy binary flags still work --------------------------------------

def test_v001_legacy_flags_unchanged(v001_mask):
    # Mixed-case flag names (Cloud flag / Dilated Cloud Flag / Cirrus flag)
    # must all be recognized as flags.
    bands = [0, 1, 2]
    qmask = quality_mask(v001_mask["path"], bands)
    _assert_binary_mask(qmask)
    np.testing.assert_array_equal(qmask, _expected_or(v001_mask["array"], bands))


def test_v001_single_aggregate_flag(v001_mask):
    qmask = quality_mask(v001_mask["path"], [7])
    _assert_binary_mask(qmask)
    np.testing.assert_array_equal(qmask, _expected_or(v001_mask["array"], [7]))


@pytest.mark.parametrize("band, name", [(5, "AOD550"), (6, "H2O (g cm-2)")])
def test_v001_data_bands_rejected(v001_mask, band, name):
    with pytest.raises(ValueError) as exc:
        quality_mask(v001_mask["path"], [band])
    assert name in str(exc.value)


def test_v001_mixed_flag_and_data_rejected(v001_mask):
    with pytest.raises(ValueError) as exc:
        quality_mask(v001_mask["path"], [0, 5])
    assert "AOD550" in str(exc.value)


# --- V002: binary flags work ---------------------------------------------------

def test_v002_binary_flags_combined(v002_mask):
    # Includes the new V002 binary SpecTf-Cloud Flag at index 9.
    bands = [0, 1, 4, 7, 9]
    qmask = quality_mask(v002_mask["path"], bands)
    _assert_binary_mask(qmask)
    np.testing.assert_array_equal(qmask, _expected_or(v002_mask["array"], bands))


def test_v002_all_zero_flag(v002_mask):
    # Water Flag (idx 2) is all zeros in this fixture.
    qmask = quality_mask(v002_mask["path"], [2])
    _assert_binary_mask(qmask)
    assert qmask.sum() == 0


# --- V002: the silent-failure case now raises ----------------------------------

@pytest.mark.parametrize(
    "band, name",
    [(8, "SpecTf-Cloud Probability"), (10, "SpecTf-Buffer Distance")],
)
def test_v002_continuous_band_rejected(v002_mask, band, name):
    with pytest.raises(ValueError) as exc:
        quality_mask(v002_mask["path"], [band])
    assert name in str(exc.value)


def test_v002_old_behavior_was_silently_wrong(v002_mask):
    """Document the bug: the previous sum-then-clip logic silently accepted a
    continuous band and returned a non-binary array, while the fixed function
    refuses it."""
    arr = v002_mask["array"]
    # Reproduce the previous implementation on the SpecTf-Cloud Probability band.
    legacy = np.sum(arr[:, :, [8]], axis=-1)
    legacy[legacy > 1] = 1
    legacy_values = set(float(v) for v in np.unique(legacy))
    # The "mask" the old code returned contained fractional probabilities, i.e.
    # it was never a valid {0, 1} mask.
    assert not legacy_values.issubset({0.0, 1.0})
    # The fixed function refuses the same band instead.
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], [8])


# --- V002: probability accepted only via explicit threshold --------------------

def test_v002_probability_threshold(v002_mask):
    qmask = quality_mask(v002_mask["path"], [8], threshold=0.5)
    _assert_binary_mask(qmask)
    expected = (v002_mask["array"][:, :, 8] >= 0.5).astype(np.uint8)
    np.testing.assert_array_equal(qmask, expected)


def test_v002_probability_threshold_combined_with_flag(v002_mask):
    qmask = quality_mask(v002_mask["path"], [0, 8], threshold=0.5)
    _assert_binary_mask(qmask)
    expected = (v002_mask["array"][:, :, 0] > 0).astype(np.uint8)
    expected |= (v002_mask["array"][:, :, 8] >= 0.5).astype(np.uint8)
    np.testing.assert_array_equal(qmask, expected)


def test_threshold_only_applies_to_probability(v002_mask):
    # Buffer Distance is continuous but not a probability: threshold must not
    # silently accept it.
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], [10], threshold=0.5)


@pytest.mark.parametrize("bad", [-0.1, 1.5])
def test_threshold_out_of_range(v002_mask, bad):
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], [8], threshold=bad)


# --- Argument validation -------------------------------------------------------

def test_out_of_range_index_rejected(v002_mask):
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], [11])  # only 0..10 exist


def test_single_int_accepted(v002_mask):
    qmask = quality_mask(v002_mask["path"], 0)
    _assert_binary_mask(qmask)
    np.testing.assert_array_equal(qmask, _expected_or(v002_mask["array"], [0]))


# --- Optional: verify against a real downloaded granule ------------------------

@pytest.mark.skipif(
    not os.environ.get("EMIT_L2A_MASK_V002"),
    reason="set EMIT_L2A_MASK_V002 to a real EMIT L2A MASK V002 granule to run",
)
def test_real_v002_granule():
    fp = os.environ["EMIT_L2A_MASK_V002"]
    # Binary flags combine into a valid mask.
    qmask = quality_mask(fp, [0, 1, 4])
    _assert_binary_mask(qmask)
    # SpecTf-Cloud Probability (index 8) is rejected unless a threshold is given.
    with pytest.raises(ValueError):
        quality_mask(fp, [8])
    qmask_thr = quality_mask(fp, [8], threshold=0.5)
    _assert_binary_mask(qmask_thr)
