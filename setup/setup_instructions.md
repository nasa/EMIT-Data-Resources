# Repository Setup Instructions

The how-tos and tutorials in this repository require a [NASA Earthdata account](https://urs.earthdata.nasa.gov/), an installation of [Git](https://git-scm.com/downloads), and a compatible Python Environment. Resources in the EMIT-Data-Resources repository have been developed using the Openscapes 2i2c JupyterHub cloud workspace.

For local Python environment setup we recommend using [mamba](https://mamba.readthedocs.io/en/latest/) to manage Python packages. To install *mamba*, download [miniforge](https://github.com/conda-forge/miniforge) for your operating system.  If using Windows, be sure to check the box to "Add mamba to my PATH environment variable" to enable use of mamba directly from your command line interface. **Note that this may cause an issue if you have an existing mamba install through Anaconda.**  

## Python Environment Setup

These Python Environments will work for all of the guides, how-to's, and tutorials within this repository, and the [VITALS repository](https://github.com/nasa/VITALS).

1. Using your preferred command line interface (command prompt, terminal, cmder, etc.) navigate to your local copy of the repository, then type the following to create a compatible Python environment.

    For Windows:

    ```cmd
    mamba create -n lpdaac_vitals -c conda-forge --yes python=3.10 fiona=1.8.22 gdal hvplot geoviews rioxarray rasterio jupyter geopandas earthaccess jupyter_bokeh h5py h5netcdf spectral scikit-image jupyterlab seaborn dask ray-default
    ```

    For MacOSX:

    ```cmd
    mamba create -n lpdaac_vitals -c conda-forge --yes python=3.10 gdal=3.7.2 hvplot geoviews rioxarray rasterio geopandas fiona=1.9.4 jupyter earthaccess jupyter_bokeh h5py h5netcdf spectral scikit-image seaborn jupyterlab dask ray-default
    ```

2. Next, activate the Python Environment that you just created.

    ```cmd
    mamba activate lpdaac_vitals 
    ```

3. Now you can launch Jupyter Notebook to open the notebooks included.

    ```cmd
    jupyter notebook 
    ```

**Still having trouble getting a compatible Python environment set up? Contact [LP DAAC User Services](https://lpdaac.usgs.gov/lpdaac-contact-us/).**  

---

## Contact Info  

Email: <LPDAAC@usgs.gov>  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 03-13-2024  

¹Work performed under USGS contract G15PD00467 for NASA contract NNG14HH33I.  
