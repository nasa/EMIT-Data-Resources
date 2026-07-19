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
import pytest
import xarray as xr

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


def _write_mask_nc(path, band_names, mask_array, fill_value=None):
    """Write a minimal EMIT-like L2A Mask netCDF file with the h5netcdf backend.

    Structure mirrors the real product: a root ``mask`` variable with dims
    ``(downtrack, crosstrack, bands)`` and a ``sensor_band_parameters`` group
    holding the ``mask_bands`` string variable. Uses only xarray + h5netcdf (both
    already required by emit_tools), so the tests add no netCDF backend
    dependency. When ``fill_value`` is given, the ``mask`` variable declares that
    ``_FillValue`` (the real V002 product declares ``_FillValue = -9999``, which
    xarray decodes to NaN on read).
    """
    path = str(path)
    root = xr.Dataset(
        {"mask": (("downtrack", "crosstrack", "bands"), mask_array.astype("float32"))}
    )
    root["mask"].attrs = {"long_name": "Masks", "units": "unitless"}
    encoding = None if fill_value is None else {"mask": {"_FillValue": fill_value}}
    root.to_netcdf(path, engine="h5netcdf", encoding=encoding)
    sbp = xr.Dataset({"mask_bands": (("bands",), np.array(band_names))})
    sbp["mask_bands"].attrs = {"long_name": "Mask Band Names"}
    sbp.to_netcdf(path, engine="h5netcdf", group="sensor_band_parameters", mode="a")
    return path


class _DatasetSpy:
    """Transparent proxy around an xarray Dataset that records ``close()``.

    Used to prove ``quality_mask`` closes the file handles it opens (including
    on error) without depending on private xarray internals.
    """

    def __init__(self, ds):
        self._ds = ds
        self.closed = False

    def __getattr__(self, item):
        return getattr(self._ds, item)

    def __getitem__(self, item):
        return self._ds[item]

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def close(self):
        self.closed = True
        return self._ds.close()


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


@pytest.mark.parametrize("bad", [True, False, np.bool_(True), np.bool_(False)])
def test_bool_index_rejected(v002_mask, bad):
    # bool is a subclass of int; it must not be accepted as a band index.
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], [bad])


def test_bool_scalar_index_rejected(v002_mask):
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], True)


def test_axis_metadata_mismatch_raises(tmp_path):
    # 'mask' has 11 bands but 'mask_bands' lists only 8 -> inconsistent layout.
    path = str(tmp_path / "mismatch.nc")
    root = xr.Dataset(
        {"mask": (("downtrack", "crosstrack", "bands"), np.zeros((3, 4, 11), "float32"))}
    )
    root.to_netcdf(path, engine="h5netcdf")
    sbp = xr.Dataset({"mask_bands": (("mask_band",), np.array(V001_MASK_BANDS))})
    sbp.to_netcdf(path, engine="h5netcdf", group="sensor_band_parameters", mode="a")
    with pytest.raises(ValueError):
        quality_mask(path, [0])


# --- Non-finite (fill / no-data) handling --------------------------------------

def test_flag_nan_pixel_is_masked_not_clear(tmp_path):
    # A no-data (NaN) pixel in a flag layer must be EXCLUDED (1), never returned
    # as a clean observation (0). _A[2, 3] is 0 (clear) before we blank it.
    cloud = _A.copy()
    cloud[2, 3] = np.nan
    layers = [cloud, _B, _ZERO, _ZERO, _C, _AOD, _H2O, _D, _PROB, _E, _DIST]
    arr = np.stack(layers, axis=-1)
    path = _write_mask_nc(tmp_path / "v002_nan.nc", V002_MASK_BANDS, arr)

    qmask = quality_mask(path, [0])
    _assert_binary_mask(qmask)
    assert qmask[2, 3] == 1
    # The previous `(layer > 0)` logic would have returned 0 (clear) here.
    assert bool((np.nan_to_num(arr[:, :, 0]) > 0)[2, 3]) is False


def test_flag_fillvalue_pixel_is_masked(tmp_path):
    # Same, but via a declared _FillValue (-9999) as in the real product, which
    # xarray decodes to NaN on read. _A[1, 0] is 0 (clear) before we blank it.
    cloud = _A.copy()
    cloud[1, 0] = -9999.0
    layers = [cloud, _B, _ZERO, _ZERO, _C, _AOD, _H2O, _D, _PROB, _E, _DIST]
    arr = np.stack(layers, axis=-1)
    path = _write_mask_nc(
        tmp_path / "v002_fill.nc", V002_MASK_BANDS, arr, fill_value=-9999.0
    )

    qmask = quality_mask(path, [0])
    _assert_binary_mask(qmask)
    assert qmask[1, 0] == 1


def test_flag_infinite_pixels_are_masked(tmp_path):
    # Both +Inf and -Inf are non-finite and must be excluded (1). Naive
    # `layer > 0` would return -Inf as 0 (clear) -- the fail-closed bug.
    # _A[0, 1] and _A[2, 0] are both 0 (clear) before we blank them.
    cloud = _A.copy()
    cloud[0, 1] = np.inf
    cloud[2, 0] = -np.inf
    layers = [cloud, _B, _ZERO, _ZERO, _C, _AOD, _H2O, _D, _PROB, _E, _DIST]
    arr = np.stack(layers, axis=-1)
    path = _write_mask_nc(tmp_path / "v002_inf.nc", V002_MASK_BANDS, arr)

    qmask = quality_mask(path, [0])
    _assert_binary_mask(qmask)
    assert qmask[0, 1] == 1  # +inf excluded
    assert qmask[2, 0] == 1  # -inf excluded


def test_threshold_nan_pixel_is_masked(tmp_path):
    # A NaN in the probability layer must also be excluded, not clear.
    # _PROB[0, 0] is 0.30, which is < 0.5 (clear) before we blank it.
    prob = _PROB.copy()
    prob[0, 0] = np.nan
    layers = [_A, _B, _ZERO, _ZERO, _C, _AOD, _H2O, _D, prob, _E, _DIST]
    arr = np.stack(layers, axis=-1)
    path = _write_mask_nc(tmp_path / "v002_prob_nan.nc", V002_MASK_BANDS, arr)

    qmask = quality_mask(path, [8], threshold=0.5)
    _assert_binary_mask(qmask)
    assert qmask[0, 0] == 1
    # `(layer >= 0.5)` alone would have returned 0 (clear) here.
    assert bool((np.nan_to_num(arr[:, :, 8]) >= 0.5)[0, 0]) is False


# --- Resource management: opened datasets are closed ---------------------------

def _install_open_spy(monkeypatch):
    import emit_tools

    real_open = emit_tools.xr.open_dataset
    spies = []

    def spy_open(*args, **kwargs):
        spy = _DatasetSpy(real_open(*args, **kwargs))
        spies.append(spy)
        return spy

    monkeypatch.setattr(emit_tools.xr, "open_dataset", spy_open)
    return spies


def test_datasets_closed_on_success(v002_mask, monkeypatch):
    spies = _install_open_spy(monkeypatch)
    quality_mask(v002_mask["path"], [0, 1])
    assert len(spies) == 2
    assert all(s.closed for s in spies)


def test_datasets_closed_on_error(v002_mask, monkeypatch):
    spies = _install_open_spy(monkeypatch)
    with pytest.raises(ValueError):
        quality_mask(v002_mask["path"], [8])  # continuous band -> raises inside `with`
    assert len(spies) == 2
    assert all(s.closed for s in spies)


# --- Optional: verify against a real downloaded granule ------------------------

@pytest.mark.skipif(
    not os.environ.get("EMIT_L2A_MASK_V002"),
    reason="set EMIT_L2A_MASK_V002 to a real EMIT L2A MASK V002 granule to run",
)
def test_real_v002_granule():
    fp = os.environ["EMIT_L2A_MASK_V002"]

    # The verified V002 SpecTf layers sit at indices 8, 9, 10.
    with xr.open_dataset(fp, engine="h5netcdf", group="sensor_band_parameters") as sbp:
        names = [str(x) for x in np.asarray(sbp["mask_bands"].data).ravel()]
    assert len(names) == 11
    assert names[8] == "SpecTf-Cloud Probability"
    assert names[9] == "SpecTf-Cloud Flag"
    assert names[10] == "SpecTf-Buffer Distance"

    # Binary flags combine into a valid mask.
    qmask = quality_mask(fp, [0, 1, 4])
    _assert_binary_mask(qmask)

    # SpecTf-Cloud Probability (index 8) is rejected unless a threshold is given.
    with pytest.raises(ValueError):
        quality_mask(fp, [8])

    # With a threshold, the result equals a directly-computed fail-closed
    # reference over the decoded probability layer.
    t = 0.5
    with xr.open_dataset(fp, engine="h5netcdf") as ds:
        raw8 = np.asarray(ds["mask"].isel(bands=8).values)
    reference = np.where(np.isfinite(raw8), raw8 >= t, True).astype(np.uint8)
    qmask_thr = quality_mask(fp, [8], threshold=t)
    _assert_binary_mask(qmask_thr)
    np.testing.assert_array_equal(qmask_thr, reference)
