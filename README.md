# EMIT-Data-Resources  

<img alt="Python resources and guides for working with EMIT data." src="img/earthdata_logo.png" width="50%" align="center" />

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![Cite this repo](https://img.shields.io/badge/Cite-EMIT--Data--Resources-blue)](https://github.com/nasa/EMIT-Data-Resources/blob/main/CITATION.cff)

Welcome to the EMIT-Data-Resources repository. This repository provides guides, short how-tos, and tutorials to help users access and work with data from the [Earth Surface Mineral Dust Source Investigation (EMIT) mission](https://lpdaac.usgs.gov/data/get-started-data/collection-overview/missions/emit-overview/). In the interest of open science this repository has been made public but is still under active development. All notebooks and scripts should be functional, however, changes or additions may be made. Make sure to consult the [CHANGE_LOG.md](CHANGE_LOG.md) for the most recent changes to the repository. Contributions from all parties are welcome.  

---

## EMIT Background  

The [EMIT](https://earth.jpl.nasa.gov/emit/) Project delivers space-based measurements of surface mineralogy of the Earth’s arid dust source regions. These measurements are used to initialize the compositional makeup of dust sources in Earth System Models (ESMs). The dust cycle, which describe the generation, lofting, transport, and deposition of mineral dust, plays an important role in ESMs. Dust composition is presently the largest uncertainty factor in quantifying the magnitude of aerosol direct radiative forcing. By understanding the composition of mineral dust sources, EMIT aims to constrain the sign and magnitude of dust-related radiative forcing at regional and global scales. During its one-year mission on the International Space Station (ISS), EMIT will make measurements over the sunlit Earth’s dust source regions that fall within ±52° latitude. EMIT will schedule up to five visits (three on average) of each arid target region and only acquisitions not dominated by cloud cover will be downlinked. EMIT-based maps of the relative abundance of source minerals will advance the understanding of the current and future impacts of mineral dust in the Earth system.  

EMIT Data Products are distributed by the [LP DAAC](https://www.earthdata.nasa.gov/centers/lp-daac). Learn more about EMIT data products from [EMIT Product Pages](https://www.earthdata.nasa.gov/data/catalog?keyword=%22EMIT%22) and search for and download EMIT data products using [NASA EarthData Search](https://search.earthdata.nasa.gov/search?q=%22EMIT%22)  

---

## Prerequisites/Setup Instructions  

This repository requires that users set up a compatible Python environment and download the EMIT granules used. See the `setup_instuctions.md` file in the `./setup/` folder.  

## Repository Contents  

Below are the resources available for EMIT Data.  

|Name|Type|Summary|
|:---|:---|:---|
|[Getting EMIT Data using EarthData Search](guides/Getting_EMIT_Data_using_EarthData_Search.md)|Markdown Guide|A thorough walkthrough for using [EarthData Search](https://search.earthdata.nasa.gov/search) to find and download EMIT data|
|[Streaming NASA Earthdata Cloud-Optimized GeoTIFFs using QGIS](guides/Streaming_cloud_optimized_geotiffs_using_QGIS.md)|Markdown Guide|A walkthrough to set up QGIS to stream cloud-optimized geotiff files from NASA Earthdata|
|[Exploring EMIT L2A Reflectance](python/tutorials/Exploring_EMIT_L2A_Reflectance.ipynb)|Jupyter Notebook|Explore EMIT L2A Reflectance data using interactive plots|
|[Visualizing Methane Plume Timeseries](python/tutorials/Visualizing_Methane_Plume_Timeseries.ipynb)|Jupyter Notebook|Find EMIT L2B CH4 Plume Data and build a timeseries of CH4 plume complexes|
|[Generating_Methane_Spectral_Fingerprint](python/tutorials/Generating_Methane_Spectral_Fingerprint.ipynb)|Jupyter Notebook|Extract Radiance Spectra and build an in-plume/out-of-plume ratio to compare with CH4 absorption coefficient|
|[Finding_EMIT_L2B_Mineral_Data](python/tutorials/Finding_EMIT_L2B_Mineral_Data.ipynb)|Jupyter Notebook|Use the `earthaccess` Python library to find EMIT L2B Mineral Identification Band Depth and Uncertainty data|
|[Working with EMIT L2B Mineralogy](python/tutorials/Working_with_EMIT_L2B_Mineralogy.ipynb)|Jupyter Notebook|Work with the EMIT L2B Mineral Identification Band Depth and Uncertainty Data and aggregate individual spectral library constituents into the EMIT-10 minerals and estimate abundance| 
|[How to find and access EMIT data](python/how-tos/How_to_find_and_access_EMIT_data.ipynb)|Jupyter Notebook|Use the `earthaccess` Python library to find and download or stream EMIT data|
|[How to Convert to ENVI Format](python/how-tos/How_to_Convert_to_ENVI.ipynb)|Jupyter Notebook|Convert from downloaded netCDF4 (.nc) format to .envi format|
|[How to Orthorectify](python/how-tos/How_to_Orthorectify.ipynb)|Jupyter Notebook|Use the geometry lookup table (GLT) included with the EMIT netCDF4 file to project on a geospatial grid (EPSG:4326)|
|[How to Extract Point Data](python/how-tos/How_to_Extract_Points.ipynb)|Jupyter Notebook|Extract spectra using lat/lon coordinates from a .csv and build a dataframe/.csv output|
|[How to Extract Area Data](python/how-tos/How_to_Extract_Area.ipynb)|Jupyter Notebook|Extract an area defined by a .geojson or shapefile|
|[How to use EMIT Quality Data](python/how-tos/How_to_use_EMIT_Quality_data.ipynb)|Jupyter Notebook|Build a mask using an EMIT L2A Mask file and apply it to an L2A Reflectance file|
|[How to use Direct S3 Access with EMIT](python/how-tos/How_to_Direct_S3_Access.ipynb)|Jupyter Notebook|Use S3 from inside AWS us-west2 to access EMIT Data|
|[How to find EMIT Data using NASA's CMR API](python/how-tos/How_to_find_EMIT_data_using_CMR_API.ipynb)|Jupyter Notebook|Use NASA's CMR API to programmatically find EMIT Data|

## Citation
This repository is the product of a collaborative effort from the NASA Land Processes Distributed Active Archive Center (DAAC) and NASA Jet Propulsion Laboratory. If you use this resource or an associated data product for your research, we appreciate a citation.

[![Cite this repo](https://img.shields.io/badge/Cite-EMIT--Data--Resources-blue)](https://github.com/nasa/EMIT-Data-Resources/blob/main/CITATION.cff)

For each data product used, you can retrieve a citation in your desired format from the product pages below.

|Product|Product Page|Collection Shortname|
|:---|:---|:---|
|EMIT L1B At-Sensor Calibrated Radiance and Geolocation Data|https://doi.org/10.5067/EMIT/EMITL1BRAD.001|EMITL1BRAD|
|EMIT L2A Estimated Surface Reflectance and Uncertainty and Masks|https://doi.org/10.5067/EMIT/EMITL2ARFL.001|EMITL2ARFL|
|EMIT L2B Estimated Mineral Identification and Band Depth and Uncertainty|https://doi.org/10.5067/EMIT/EMITL2BMIN.001|EMITL2BMIN|
|EMIT L2B Estimated Methane Enhancement Data|https://doi.org/10.5067/EMIT/EMITL2BCH4ENH.001|EMITL2BCH4ENH|
|EMIT L2B Estimated Methane Plume Complexes|https://doi.org/10.5067/EMIT/EMITL2BCH4PLM.001|EMITL2BCH4PLM|


## Contributors

If you would like to contribute, please view our [contributing guide](CONTRIBUTING.md).

[![Contributors](https://contrib.rocks/image?repo=nasa/emit-data-resources)](https://github.com/nasa/emit-data-resources/graphs/contributors)

## Related Resources

LP DAAC also develops and maintains additional resources, available at <https://github.com/nasa/LPDAAC-Data-Resources>. This landing repository provides general guides and tutorials, as well as links to mission-specific repositories that offer detailed guides, tutorials, how-tos, and scripts to help users find, access, and work with each mission’s data products.

| GitHub Repository | Summary | Contents Description |
|----|-----|----|
|[AppEEARS Data Resources](https://github.com/nasa/AppEEARS-Data-Resources) |How to use the Application for Extracting and Exploring Analysis Ready Samples (AppEEARS) |Tutorials, AppEEARS API, Direct S3 Access |
|[AppEEARS QGIS Plugin](https://github.com/nasa/AppEEARS-QGIS-Plugin)|Browse and load cloud-optimized geotiff output files from AppEEARS directly into QGIS|QGIS Plugin|
|[ASTER Data Resources](https://github.com/nasa/ASTER-Data-Resources)| How to find, access, and work with ASTER data |Tutorials|
|[ECOSTRESS Data Resources](https://github.com/nasa/ECOSTRESS-Data-Resources)|How to find, access, and work with ECOSTRESS data (The ECOsystem Spaceborne Thermal Radiometer Experiment on Space Station)|Tutorials, Scripts, Direct S3 Access|
|[EMIT Data Resources](https://github.com/nasa/EMIT-Data-Resources) |How to find, access, and work with EMIT data (Earth Surface Mineral Dust Source Investigation)|Tutorials, Scripts, Direct S3 Access |
|[GEDI Data Resources](https://github.com/nasa/GEDI-Data-Resources) |How to find, access, and work with GEDI data (Global Ecosystem Dynamics Investigation)|Tutorials |
|[HLS Data Resources](https://github.com/nasa/HLS-Data-Resources)|How to find, access, and work with HLS data (Harmonized Landsat Sentinel-2)|Tutorials, Scripts, Direct S3 Access|
|[LPDAAC Data Resources](https://github.com/nasa/LPDAAC-Data-Resources)|How to find, access, and work with LPDAAC data |Tutorials, Scripts, Direct S3 Access|
|[MODIS-VIIRS Data Resources](https://github.com/nasa/MODIS-VIIRS-Data-Resources)|How to find, access, and work with MODIS and VIIRS data|Tutorials|
|[VITALS](https://github.com/nasa/VITALS)|How to find and work with EMIT and ECOSTRESS data together |Tutorials|

### Other Helpful Links   

+ [NASA Earthdata Website](https://www.earthdata.nasa.gov/)
+ [LP DAAC Website](https://www.earthdata.nasa.gov/centers/lp-daac)
+ [JPL EMIT Website](https://earth.jpl.nasa.gov/emit/)  
+ [EMIT Tutorial Webinar Series Recordings](https://www.youtube.com/playlist?list=PLO2yB4LGNlWrC5NdxeHMxyAxdwQhSypXe)
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
Date last modified: 06-28-2024  

¹Work performed under USGS contract G15PD00467 for NASA contract NNG14HH33I.  
