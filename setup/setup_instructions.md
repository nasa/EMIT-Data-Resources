# Repository Setup Instructions

The how-tos and tutorials in this repository require a [NASA Earthdata account](https://urs.earthdata.nasa.gov/), an installation of [Git](https://git-scm.com/downloads), and a compatible Python Environment. To manage Python Environments we recommend using [Anaconda](https://www.anaconda.com/products/distribution) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) as a package manager. If using Windows, be sure to check the box to "Add Anaconda to my PATH environment variable" to enable use of conda directly from your command line interface. Note that this will cause an issue if you have an existing version of conda or miniconda installed.  

## Python Environment Setup

This Python Environment will work for all of the guides, how-to's, and tutorials within this repository. It may take a while to download and configure all of the packages properly.  

1. Using your preferred command line interface (command prompt, terminal, cmder, etc.) navigate to your local copy of the repository, then type the following to create a compatible Python environment.  

    ```cmd
    conda create -n emit_tutorials -c conda-forge --yes python=3.10 gdal=3.7.0 hvplot=0.8.4 geoviews=1.9.6 rioxarray rasterio jupyter geopandas earthaccess jupyter_bokeh h5py h5netcdf spectral
    ```

2. Next, activate the Python Environment that you just created.

    ```cmd
    conda activate emit_tutorials 
    ```

3. Now you can launch Jupyter Notebook to open the notebooks included.

    ```cmd
    jupyter notebook 
    ```

If you're having trouble creating a compatible Python Environment, you can also create one from the included `.yml` file. Using your preferred command line interface (command prompt, terminal, cmder, etc.) navigate to your local copy of the repository, then type the following to create a compatible Python environment using the `.yml` file.

```cmd
conda env create -f setup/emit_tutorials.yml
```

After this, follow steps 2 and 3 above.

[Additional information](https://conda.io/docs/user-guide/tasks/manage-environments.html) on setting up and managing Conda environments.  
**Still having trouble getting a compatible Python environment set up? Contact [LP DAAC User Services](https://lpdaac.usgs.gov/lpdaac-contact-us/).**  

---

## Contact Info  

Email: <LPDAAC@usgs.gov>  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 01-09-2023  

¹Work performed under USGS contract G15PD00467 for NASA contract NNG14HH33I.  
