# Change Log

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
_________________________________________________________________________

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
