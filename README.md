# EMIT-Data-Resources

Welcome to the EMIT-Data-Resources repository. This repository provides guides, short how-tos, and tutorials to help users access and work with data from the Earth Surface Mineral Dust Source Investigation (EMIT) mission. In the interest of open science this repository has been made public but is still under active development. All jupyter notebooks and scripts should be functional, however, changes or additions may be made. Contributions from all parties are welcome.

These guides, how-tos, and tutorials can also be found in the [EMIT Tutorials Repository](https://github.com/emit-sds/tutorials)

---

## Background  

The EMIT Project delivers space-based measurements of surface mineralogy of the Earth’s arid dust source regions. These measurements are used to initialize the compositional makeup of dust sources in Earth System Models (ESMs). The dust cycle, which describe the generation, lofting, transport, and deposition of mineral dust, plays an important role in ESMs.  Dust composition is presently the largest uncertainty factor in quantifying the magnitude of aerosol direct radiative forcing.  By understanding the composition of mineral dust sources, EMIT aims to constrain the sign and magnitude of dust-related radiative forcing at regional and global scales. During its one-year mission on the International Space Station (ISS), EMIT will make measurements over the sunlit Earth’s dust source regions that fall within ±52° latitude. EMIT will schedule up to five visits (three on average) of each arid target region and only acquisitions not dominated by cloud cover will be downlinked. EMIT-based maps of the relative abundance of source minerals will advance the understanding of the current and future impacts of mineral dust in the Earth system.  

---

## Prerequisites/Setup Instructions

This repository requires that users set up a compatible Python environment and download the EMIT granules used. See the sections below to accomplish these prerequisites.

### 1. Python Environment Setup  

This Python Environment will work for all of the guides, how-to's, and tutorials within this repository. Using your preferred command line interface (command prompt, terminal, cmder, etc.) navigate to your local copy of the repository, then type the following to create a compatible Python environment using the included `.yml` file.  

> `conda env create -f emit_tutorials.yml`  

Next, activate the Python Environment that you just created.

> `conda activate emit_tutorials`  

Now you can launch Jupyter Notebook to open the notebooks included.

> `jupyter notebook`  

[Additional information](https://conda.io/docs/user-guide/tasks/manage-environments.html) on setting up and managing Conda environments.  
**Still having trouble getting a compatible Python environment set up? Contact [LP DAAC User Services](https://lpdaac.usgs.gov/lpdaac-contact-us/).**  

### 2. File Downloads  

These granules below are used within this tutorial. Click/copy the URLs into a browser to download. Save them into the `/data/` folder within this repository. You will need a [NASA Earth Data Search](https://search.earthdata.nasa.gov/search) login.

+ L2A Reflectance Granule - <https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2ARFL.001/EMIT_L2A_RFL_001_20220903T163129_2224611_012/EMIT_L2A_RFL_001_20220903T163129_2224611_012.nc>  
+ L2A Mask Granule - <https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2ARFL.001/EMIT_L2A_RFL_001_20220903T163129_2224611_012/EMIT_L2A_MASK_001_20220903T163129_2224611_012.nc>  

You can also choose to use the DAAC Data Downloader Tool, a tool designed to automate the download process by reading a list of URLs from a text file. Open your preferred command line interface and navigate to the your local copy of this repository. To download the repository, type the following in the command line. This requires a NASA Earth Data account and will prompt you to configure a `.netrc` file if you do not have one configured.

>`git clone https://git.earthdata.nasa.gov/scm/lpdur/daac_data_download_python.git ../daac_data_download_python/`

After copying the repository you can execute the python script, providing the `-dir` argument as a directory to save to, and the `-f` argument, a URL or text file containing a list of URLs to download.  

>`python ../daac_data_download_python/DAACDataDownload.py -dir ./data/ -f ./data/emit_data_urls.txt`

---

## Repository Contents  

### **Guides**  

+ Getting EMIT Data using EarthData Search - A thourough walkthrough for using [EarthData Search](https://search.earthdata.nasa.gov/search) to find and download EMIT data.

### **How-To Notebooks**

+ How to Convert to ENVI Format - Convert from downloaded netCDF4 (.nc) format to .envi format.
+ How to Orthorectify - Use the geometry lookup table (GLT) included with the EMIT netCDF4 file to project on a geospatial grid (EPSG:4326).
+ How to Extract Point Data  - Extract spectra using lat/lon coordinates from a .csv and build a dataframe/.csv output.
+ How to Extract Area Data - Clip to/extract an area defined by a .geojson or shapefile.
+ How to use EMIT Quality Data - Build a mask using bands from  an EMIT L2A Mask file and apply it to an L2A Reflectance file.

### **Tutorial Notebooks**  

1. Exploring EMIT L2A Reflectance  
    1.1 Setup  
    1.2 Opening and Understanding File Structure  
    1.3 Plotting Spectra Basics  
    1.4 Geocorrection/Applying GLT  
    1.5 Spatial Plotting (Imagery)  
    1.6 Writing an Orthorectified netCDF4 Output
    1.7 Exploring Spectral and Spatial Plots  

---

## Helpful Links  

+ [EMIT Website](https://earth.jpl.nasa.gov/emit/)  

+ [EMIT Github Repository](https://github.com/emit-sds) - Main EMIT Repository  

+ [EMIT Utilities Github Repository](https://github.com/emit-sds/emit-utils) - General convenience utilities for working with EMIT data

+ [L2A Reflectance User Guide](https://lpdaac.usgs.gov/documents/1569/EMITL2ARFL_User_Guide_v1.pdf)  

+ [L2A Algorithm Theoretical Basis Document](https://lpdaac.usgs.gov/documents/1571/EMITL2A_ATBD_v1.pdf)  

+ [EMIT LP DAAC Product Pages](https://lpdaac.usgs.gov/product_search/?query=emit&status=Operational&view=cards&sort=title) - Learn more about available EMIT products  

+ [EMIT on Earth Data Search](https://search.earthdata.nasa.gov/search?q=%22EMIT%22) - Find EMIT Data  

---

## Contact Info:  

Email: LPDAAC@usgs.gov  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 11-21-2022  

¹Work performed under USGS contract G15PD00467 for NASA contract NNG14HH33I.  
