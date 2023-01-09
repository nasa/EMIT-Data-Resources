# Repository Setup Instructions

The how-tos and tutorials in this repository require a compatible Python Environment, an installation of [Git](https://git-scm.com/downloads), and some EMIT granules that are too large to be included in the repository. To setup the environment and download these files, follow the steps in sections 1 and 2.  

+ If you do not have an Environment Manager installed, we recommend  [Anaconda](https://www.anaconda.com/products/distribution) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). When installing, be sure to check the box to "Add Anaconda to my PATH environment variable" to enable use of conda directly from your command line interface.
+ If you do not have Git, you can download it [here](https://git-scm.com/downloads).  

## 1. Python Environment Setup  

This Python Environment will work for all of the guides, how-to's, and tutorials within this repository. Using your preferred command line interface (command prompt, terminal, cmder, etc.) navigate to your local copy of the repository, then type the following to create a compatible Python environment using the included `.yml` file.  

> `conda env create -f setup/emit_tutorials.yml`  

Next, activate the Python Environment that you just created.

> `conda activate emit_tutorials`  

Now you can launch Jupyter Notebook to open the notebooks included.

> `jupyter notebook`  

[Additional information](https://conda.io/docs/user-guide/tasks/manage-environments.html) on setting up and managing Conda environments.  
**Still having trouble getting a compatible Python environment set up? Contact [LP DAAC User Services](https://lpdaac.usgs.gov/lpdaac-contact-us/).**  

## 2. File Downloads  

These granules below are used within this tutorial. Click/copy the URLs into a browser to download. Save them into the `./data/` folder within this repository. You will need a [NASA Earth Data Search](https://search.earthdata.nasa.gov/search) login.

+ L2A Reflectance Granule - <https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2ARFL.001/EMIT_L2A_RFL_001_20220903T163129_2224611_012/EMIT_L2A_RFL_001_20220903T163129_2224611_012.nc>  
+ L2A Mask Granule - <https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2ARFL.001/EMIT_L2A_RFL_001_20220903T163129_2224611_012/EMIT_L2A_MASK_001_20220903T163129_2224611_012.nc>  

If you have a textfile of urls you wish to download, you can also choose to use the DAAC Data Downloader Tool, a tool designed to automate the download process by reading a list of URLs from a text file. Open your preferred command line interface and navigate to the your local copy of this repository. To download the repository, type the following in the command line. This requires a NASA Earth Data account and will prompt you to configure a `.netrc` file if you do not have one configured.

>`git clone https://git.earthdata.nasa.gov/scm/lpdur/daac_data_download_python.git ../daac_data_download_python/`

After copying the repository you can execute the python script, providing the `-dir` argument as a directory to save to, and the `-f` argument, a URL or text file containing a list of URLs to download.  

>`python ../daac_data_download_python/DAACDataDownload.py -dir ./data/ -f ./data/emit_data_urls.txt`

---

## Contact Info:  

Email: LPDAAC@usgs.gov  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 01-09-2023  

¹Work performed under USGS contract G15PD00467 for NASA contract NNG14HH33I.  
