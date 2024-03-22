# Change Log

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
_________________________________________________________________________

## 2024-03-22

> ### Added
>
> - `geojson` with plume bounding box for `Visualizing Methane Plume Timeseries` notebook

> ### Changed
>
> - Updates to instructions in `Visualizing Methane Plume Timeseries` notebook
> - Remove horizontal lines from notebooks for web-book creation

## 2024-03-13

> ### Added
>
> - Added new `Visualizing Methane Plume Timeseries` notebook  
> - Added new `Generating Methane Spectral Fingerprint` notebook
> - Added new `tutorial_utils.py` module which has functions specific to the new CH4 notebooks  

> ### Changed
>
> - Updated `setup_instructions.md` to provide more up to date Python setup instructions  
> - Updated `README.md` to include new notebooks  
> - Added some functions to `emit_tools.py` to support the new methane tutorial notebooks
> - Minor updates to `How to Extract Area` and `How to Extract Points` notebooks to improve visualizations

## 2023-12-01

> ### Changed
>
> - Applied automatic formatting to `emit_tools.py`  
> - Fixed GLT index adjustment in `raw_spatial_crop` function
> - Reimplemented `PointerXY` stream so interactive plot from `Exploring EMIT L2A Reflectance` shows spectra at current cursor location

## 2023-11-27

> ### Changed
>
> - Updated `emit_tools.py` to fix an issue with `raw_spatial_crop`
> - Updated some plotting visuals within notebooks
> - Improved interactive plots in `Exploring_EMIT_L2A_Reflectance`
> - Renamed `CONTRIBUTE.md` to `CONTRIBUTING.md`

## 2023-10-23

> ### Changed
>
> - Updated `emit_tools.py` to fix an issue with application of the GLT to elevation

## 2023-10-12

> ### Changed
>
> - Updated data download method in how-to and tutorial notebooks

## 2023-08-04

> ### Changed
>
> - Updated `How to Direct S3 Access` notebook to include `earthaccess`

## 2023-07-06

> ### Added
>
> - `How_to_find_and_access_EMIT_data` Notebook - A guide for how to use `earthaccess` Python library to find and open/download EMIT data.

> ### Changed
>
> - **Default behavior of `emit_xarray` function from `emit_tools.py` no longer orthorectifies, must add parameter `ortho=True`**
> - Removed 'bands' dimension in favor of 'wavelengths' for 'radiance' and 'reflectance' when using `emit_xarray` function from `emit_tools.py`
> - Updated `write_envi` function from `emit_tools.py`, it should work for L1A,L1B, and L2A products.
> - various other minor changes to `emit_tools.py`
> - Updated notebooks to support changes made to `emit_tools.py`
> - Improved RGB Image creation and brightness adjustment in `Exploring_EMIT_L2A_Reflectance` Tutorial
> - Each notebook now uses the `earthaccess` python library to handle NASA Earthdata Login and downloads of needed files
> - Updated `setup_instructions.md` and recommended python environment `emit_tutorials_windows.yml`
> - Updated `setup_instructions.md` to recommend mamba

## 2023-05-30

> ### Changed
>
> - Restructured repository layout to match LP DAAC standards

## 2023-05-16
>
> ### Changed
>
> - Updated `emit_tools.py` to work for HTTPFileSystem URIs

## 2023-04-11
>
> ### Changed
>
> - The orthorectification process within `emit_tools.py` was creating a lat/lon grid based upon the upper left pixel coordinates. Corrected this to be pixel centers.

## 2023-03-29
  
> ### Added
>
> - CHANGE_LOG.md
> - CODE_OF_CONDUCT.md
> - CONTRIBUTING.md
> - LICENSE.md
>
> ### Changed
>
> - content in CHANGE_LOG.md
>
> ### Fixed
>
> -
