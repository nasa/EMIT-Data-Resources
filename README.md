# EMIT-Data-Resources  

Welcome to the EMIT-Data-Resources repository. This repository provides guides, short how-tos, and tutorials to help users access and work with data from the [Earth Surface Mineral Dust Source Investigation (EMIT) mission](https://lpdaac.usgs.gov/data/get-started-data/collection-overview/missions/emit-overview/). In the interest of open science this repository has been made public but is still under active development. All notebooks and scripts should be functional, however, changes or additions may be made. Make sure to consult the [CHANGE_LOG.md](CHANGE_LOG.md) for the most recent changes to the repository. Contributions from all parties are welcome.  

---

## EMIT Background  

The [EMIT](https://earth.jpl.nasa.gov/emit/) Project delivers space-based measurements of surface mineralogy of the Earth’s arid dust source regions. These measurements are used to initialize the compositional makeup of dust sources in Earth System Models (ESMs). The dust cycle, which describe the generation, lofting, transport, and deposition of mineral dust, plays an important role in ESMs. Dust composition is presently the largest uncertainty factor in quantifying the magnitude of aerosol direct radiative forcing. By understanding the composition of mineral dust sources, EMIT aims to constrain the sign and magnitude of dust-related radiative forcing at regional and global scales. During its one-year mission on the International Space Station (ISS), EMIT will make measurements over the sunlit Earth’s dust source regions that fall within ±52° latitude. EMIT will schedule up to five visits (three on average) of each arid target region and only acquisitions not dominated by cloud cover will be downlinked. EMIT-based maps of the relative abundance of source minerals will advance the understanding of the current and future impacts of mineral dust in the Earth system.  

EMIT Data Products are distributed by the [LP DAAC](https://lpdaac.usgs.gov/). Learn more about EMIT data products from [EMIT Product Pages](https://lpdaac.usgs.gov/product_search/?query=emit&status=Operational&view=cards&sort=title) and search for and download EMIT data products using [NASA EarthData Search](https://search.earthdata.nasa.gov/search?q=%22EMIT%22)  

---

## Prerequisites/Setup Instructions  

This repository requires that users set up a compatible Python environment and download the EMIT granules used. See the `setup_instuctions.md` file in the `./setup/` folder.  

## Repository Contents  

Below are the resources available for EMIT Data.  

|Name|Type|Summary|
|:---|:---|:---|
|[Getting EMIT Data using EarthData Search](guides/Getting_EMIT_Data_using_EarthData_Search.md)|Markdown Guide|A thorough walkthrough for using [EarthData Search](https://search.earthdata.nasa.gov/search) to find and download EMIT data|
|[Exploring EMIT L2A Reflectance](python/tutorials/Exploring_EMIT_L2A_Reflectance.ipynb)|Jupyter Notebook|Explore EMIT L2A Reflectance data using interactive plots|
|[How to find and access EMIT data](python/how-tos/How_to_find_and_access_EMIT_data.ipynb)|Jupyter Notebook|Use the `earthaccess` Python library to find and download or stream EMIT data|
|[How to Convert to ENVI Format](python/how-tos/How_to_Convert_to_ENVI.ipynb)|Jupyter Notebook|Convert from downloaded netCDF4 (.nc) format to .envi format|
|[How to Orthorectify](python/how-tos/How_to_Orthorectify.ipynb)|Jupyter Notebook|Use the geometry lookup table (GLT) included with the EMIT netCDF4 file to project on a geospatial grid (EPSG:4326)|
|[How to Extract Point Data](python/how-tos/How_to_Extract_Points.ipynb)|Jupyter Notebook|Extract spectra using lat/lon coordinates from a .csv and build a dataframe/.csv output|
|[How to Extract Area Data](python/how-tos/How_to_Extract_Area.ipynb)|Jupyter Notebook|Extract an area defined by a .geojson or shapefile|
|[How to use EMIT Quality Data](python/how-tos/How_to_use_EMIT_Quality_data.ipynb)|Jupyter Notebook|Build a mask using an EMIT L2A Mask file and apply it to an L2A Reflectance file|
|[How to use Direct S3 Access with EMIT](python/how-tos/How_to_Direct_S3_Access.ipynb)|Jupyter Notebook|Use S3 from inside AWS us-west2 to access EMIT Data|
|[How to find EMIT Data using NASA's CMR API](python/how-tos/How_to_find_EMIT_data_using_CMR_API.ipynb)|Jupyter Notebook|Use NASA's CMR API to programmatically find EMIT Data|

---

## Helpful Links  

+ [JPL EMIT Website](https://earth.jpl.nasa.gov/emit/)  
+ [Video of 2023 Tutorial Series](https://www.youtube.com/playlist?list=PLO2yB4LGNlWrC5NdxeHMxyAxdwQhSypXe)
+ [LP DAAC EMIT Product Pages](https://lpdaac.usgs.gov/product_search/?query=emit&status=Operational&view=cards&sort=title) - Learn more about available EMIT products  
+ [VISIONS Open Data Portal](https://earth.jpl.nasa.gov/emit/data/data-portal/coverage-and-forecasts/) - Learn about current and forecasted EMIT coverage  

+ [EMIT on Earth Data Search](https://search.earthdata.nasa.gov/search?q=%22EMIT%22) - Download EMIT Data from NASA

+ [EMIT Github Repository](https://github.com/emit-sds) - Main EMIT Repository  

+ [EMIT Utilities Github Repository](https://github.com/emit-sds/emit-utils) - General convenience utilities for working with EMIT data

+ [L2A Reflectance User Guide](https://lpdaac.usgs.gov/documents/1569/EMITL2ARFL_User_Guide_v1.pdf)  

+ [L2A Algorithm Theoretical Basis Document](https://lpdaac.usgs.gov/documents/1571/EMITL2A_ATBD_v1.pdf)  

+ [EMIT on Slack]( https://forms.gle/XefLVG6e6A7ezwpY9) - Join the EMIT slack community!

---

## Contact Info  

Email: <LPDAAC@usgs.gov>  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 07-07-2023  

¹Work performed under USGS contract G15PD00467 for NASA contract NNG14HH33I.  
