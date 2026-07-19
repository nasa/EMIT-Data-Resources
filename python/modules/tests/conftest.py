"""Test configuration for the ``emit_tools`` module tests.

``emit_tools`` imports the full geospatial stack (``osgeo.gdal``, ``geopandas``,
``rasterio``, ``rioxarray``, ``scikit-image``, ``spectral`` ...) at module load
time. None of that is needed to exercise ``quality_mask``, which only uses
``numpy`` and ``xarray``. To let these unit tests run in a minimal environment
(e.g. CI without GDAL installed), any of those optional dependencies that are
not importable are replaced with lightweight stand-ins before ``emit_tools`` is
imported. When the full environment is present, the real modules are used and
nothing is stubbed.
"""

import importlib
import sys
from pathlib import Path
from unittest import mock

# Make ``import emit_tools`` resolve, matching the notebooks' convention of
# adding ``python/modules`` to ``sys.path``.
MODULES_DIR = Path(__file__).resolve().parents[1]
if str(MODULES_DIR) not in sys.path:
    sys.path.insert(0, str(MODULES_DIR))

# Optional heavy dependencies that ``emit_tools`` imports but ``quality_mask``
# does not need. Stub only the ones that are actually missing.
_OPTIONAL_DEPENDENCIES = (
    "osgeo",
    "spectral",
    "spectral.io",
    "skimage",
    "geopandas",
    "rasterio",
    "rioxarray",
    "rioxarray.merge",
    "s3fs",
)

for _name in _OPTIONAL_DEPENDENCIES:
    try:
        importlib.import_module(_name)
    except ModuleNotFoundError as _exc:
        # Only stub the specific optional dependency that is genuinely absent.
        # A ModuleNotFoundError naming an unrelated module (e.g. a broken
        # transitive dependency or an ABI/version mismatch) is re-raised so it is
        # not silently masked as a passing test run.
        _missing = _exc.name or ""
        if _missing == _name or _name.startswith(_missing + "."):
            sys.modules[_name] = mock.MagicMock()
        else:
            raise
