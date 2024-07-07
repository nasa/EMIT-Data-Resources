# Streaming NASA Earthdata Cloud-Optimized Geotiffs using QGIS

**Note:** Users must have [NASA Earthdata Login](https://urs.earthdata.nasa.gov/home) to stream data.

QGIS can be used to stream NASA Earthdata cloud-optimized geotiff files. To do this QGIS uses https and the vsicurl virtual file system handler. This guide will show you how to configure QGIS to do this. There are 4 steps needed.

1. Create a `.netrc` file
2. Add Custom Environment Variables
3. Restart QGIS
4. Add Data

## 1. Create a .netrc file

Under the hood QGIS uses gdal and vsicurl to stream data via https. For NASA Earthdata, this requires a `.netrc` file to store your Earthdata login credentials. If you already have one, you can skip this step. This file should be located in your home directory. Choose one of the methods below to create the file and enter your credentials.

### A. Manual Set Up

- Download the [.netrc template file](https://github.com/nasa/LPDAAC-Data-Resources/tree/main/data/.netrc) and save it in your *home/user/* directory, where *user* is your personal user directory. For example: `C:\Users\user\.netrc` or `home/user/.netrc.`  
- Open the `.netrc` file in a text editor and replace <USERNAME> with your NASA Earthdata Login username and <PASSWORD> with your NASA Earthdata Login password.

After editing, the file should look something like this:

![Example .netrc 1](../img/example_netrc1.png)

or you can also have everything on a single line separated by spaces, like:

![example .netrc 2](../img/example_netrc2.png)

### B. Command Line  

**For Linux/MacOS:**

To Create a .netrc file, enter the following in the command line, replacing <USERNAME> and <PASSWORD> with your NASA Earthdata username and password. This will create a file in your home directory or append your NASA credentials to an existing file.

```bash
echo "machine urs.earthdata.nasa.gov login <USERNAME> password <PASSWORD>" >>~/.netrc
```

**For Windows:**

To Create a .netrc file, enter the following in the command line, replacing <USERNAME> and <PASSWORD> with your NASA Earthdata username and password. This will create a file in your home directory or append your NASA credentials to an existing file.

```cmd
echo machine urs.earthdata.nasa.gov login <USERNAME> password <PASSWORD> >> %userprofile%\.netrc
```

You can verify that the file is correct by opening with a text editor. It should look like an example in one of the figures above.

## 2. Environment Settings

Next, while still in the **Settings** menu, select **System** on the left-hand side, then scroll down to **Environment**. Select the check box and enter the Variables and Values in the table below by clicking the plus sign to add a new variable. These variables, set up a place for cookies to be stored and read from, prevent GDAL from reading all files in the directory, specify the extensions GDAL is allowed to access over HTTPS using the vsicurl virtual file system handler, and allow use of unsafe SSL connections (use with caution).

|Variable|Value|
| ----------- | ----------- |
|GDAL_HTTP_COOKIEFILE|~/cookies.txt|
|GDAL_HTTP_COOKIEJAR|~/cookies.txt|
|GDAL_DISABLE_READDIR_ON_OPEN|EMPTY_DIR|
|CPL_VSIL_CURL_ALLOWED_EXTENSIONS|TIF|
|GDAL_HTTP_UNSAFESSL|YES|

> **Note: These settings may affect streaming from other data sources.**

![Environment Settings](../img/environment_settings.png)

## 3. Restart QGIS

Restart QGIS to load the new environment settings.

## 4. Add Some Data

To add some raster data, select the add raster data button from the toolbar, select 'Protocol: HTTP(S), cloud, etc' as the source type, then enter a URI for a cloud-optimized geotiff file and press the add button. For example:

`https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/EMITL2BCH4ENH.001/EMIT_L2B_CH4ENH_001_20240423T074559_2411405_018/EMIT_L2B_CH4ENH_001_20240423T074559_2411405_018.tif`

![Add raster data](../img/add_data.png)

This should add the data to your map, and will work with any NASA Earthdata hosted cloud-optimized geotiff file.

![Added Scene](../img/example_scene.png)

>**If these instructions do not work, please verify that your `.netrc` file has the correct username and password, and is formatted as shown in [Section 1](#1-create-a-netrc-file).**

## Contact Info  

Email: <LPDAAC@usgs.gov>  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 06-25-2024  

¹Work performed under USGS contract 140G0121D0001 for NASA contract NNG14HH33I.  
