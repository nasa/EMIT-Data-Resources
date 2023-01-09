"""
This Module has the functions related to working with an EMIT dataset. This includes doing things
like opening and flattening the data to work in xarray, orthorectification, and extracting point and area samples.

Author: Erik Bolch, ebolch@contractor.usgs.gov 

Last Updated: 01/04/2022

TO DO: 
- Include a way to work with band_mask packed bands from L2A Mask
- Rework inclusion of elevation as separate function?
- Clear out redundant code and streamline
- Investigate reducing memory usage
- Ensure works with s3 access (s3fs context manager doesn't support netCDF4 package)
"""
# Packages used
import netCDF4 as nc
from osgeo import gdal
import numpy as np
import math
import pandas as pd
import xarray as xr
import rasterio as rio

def emit_xarray(filepath, ortho=True, quality_mask=None): #quality_mask=False
    """
        This function utilizes other functions in this module to streamline opening an EMIT dataset as an xarray.Dataset.
        
        Parameters:
        filepath: a filepath to an EMIT netCDF file
        ortho: True or False, whether to orthorectify the dataset or leave in crosstrack/downtrack coordinates.
        quality_mask: an array output from the build_mask function used to mask pixels based on quality flags selected in that function. 
                        Any non-orthorectified array with the proper crosstrack and downtrack dimensions can also be used.

        Returns:
        xr_ds: an xarray.Dataset constructed based on the parameters provided

        """
    # Read in Data as Xarray Datasets
    ds = xr.open_dataset(filepath,engine = 'h5netcdf')
    loc = xr.open_dataset(filepath, engine = 'h5netcdf', group='location')
    wvl = xr.open_dataset(filepath, engine = 'h5netcdf', group='sensor_band_parameters')
        
    if ortho is True:
        # Apply Orthorectification Function
        xr_ds = ortho_x(ds,loc,wvl,quality_mask)
    else:
        # Building Flat Dataset from Components
        data_vars = {**ds.variables, **loc.variables, **wvl.variables} 
        coords = {**ds.coords}
        xr_ds = xr.Dataset(data_vars=data_vars, coords = coords, attrs= ds.attrs)
        if quality_mask is not None:
            xr_ds['quality_mask'] = (list(ds.dims)[0:2], quality_mask)
    
    return xr_ds

# Function to Calculate the Lat and Lon Vectors/Coordinate Grid
def coord_vects(ds,loc):
    """
    This function calculates the Lat and Lon Coordinate Vectors using the GLT and Metadata from an EMIT dataset read into xarray.
    
    Parameters:
    ds: an xarray.Dataset containing the root variable and metadata of an EMIT dataset
    loc: an xarray.Dataset containing the 'location' group of an EMIT dataset

    Returns:
    lon, lat (numpy.array): longitute and latitude array grid for the dataset

    """
    # Retrieve Geotransform from Metadata
    GT = ds.geotransform
    # Create Array for Lat and Lon and fill
    dim_x = loc.glt_x.shape[1]
    dim_y = loc.glt_x.shape[0]
    lon = np.zeros(dim_x)
    lat = np.zeros(dim_y)
    for x in np.arange(dim_x):
        y=0 # No Rotation
        x_geo = GT[0] + x * GT[1] + y * GT[2]
        lon[x] = x_geo
    for y in np.arange(dim_y):
        x=0 # No Rotation
        y_geo = GT[3] + x * GT[4] + y * GT[5]
        lat[y] = y_geo
    return lon,lat

# Function to Apply the GLT to an array
def apply_glt(ds_array,glt_array,fill_value=-9999,GLT_NODATA_VALUE=0):
    """
    This function applies a numpy array of the EMIT glt to a numpy array of the desired dataset (i.e. reflectance, radiance, etc. 
    
    Parameters:
    ds_array: numpy array of the desired variable
    glt_array: a GLT array constructed from EMIT GLT data
    
    Returns: 
    out_ds: a numpy array of orthorectified data.
    """
    
    # Build Output Dataset
    out_ds = np.zeros((glt_array.shape[0], glt_array.shape[1], ds_array.shape[-1]), dtype=np.float32) + fill_value
    valid_glt = np.all(glt_array != GLT_NODATA_VALUE, axis=-1)
    # Adjust for One based Index
    glt_array[valid_glt] -= 1 
    out_ds[valid_glt, :] = ds_array[glt_array[valid_glt, 1], glt_array[valid_glt, 0], :]
    return out_ds

def ortho_x(ds,loc,wvl,quality_mask=None):
    """
    This function applies a geometry lookup table to the desired EMIT datasets using xarray datasets built from its root variables, 'location' group and 'sensor_band_parameters' group and returns an xarray dataset.

    Parameters:
    ds: xarray dataset of the root group
    loc: xarray dataset of the 'location' group
    wvl: xarray dataset of the 'sensor_band_parameters' group

    Returns: 
    out_xr: an xarray.Dataset containing the geolocated EMIT dataset and metadata
    """
    GLT_NODATA_VALUE = 0
    fill_value = -9999
    var_list = list(ds.variables)
    
    # Define GLT Dataset
    glt_ds = np.nan_to_num(np.stack([loc['glt_x'].data,loc['glt_y'].data],axis=-1),nan=GLT_NODATA_VALUE).astype(int)
    
    # Calculate Lat and Lon Vectors
    
    lon, lat = coord_vects(ds,loc) # Reorder this function to make sense in case of multiple variables
    
    # Create Coordinate Dictionary
    
    #elev_ds = apply_glt(loc['elev'].data[:,:,np.newaxis],glt_ds)
    coords = {'latitude':(['latitude'],lat), 'longitude':(['longitude'],lon), **wvl.variables}# unpack wvl to complete coordinates dictionary
    
    #elev_data_vars = {'elev':(['lat','lon'],np.squeeze(elev_ds))}      
    #out_elev = xr.Dataset(data_vars=elev_data_vars, coords=coords, attrs= ds.attrs)
    
    # Define Rawspace Dataset Variable Values (Typically Reflectance)    
    for var in var_list:
        raw_ds = ds[var].data

        # Apply Mask if Included
        if quality_mask is not None:
            raw_ds[quality_mask == 1] = fill_value
        
        # Apply GLT to dataset
        out_ds = apply_glt(raw_ds,glt_ds)
        
        #Update variables
        data_vars = {var:(['latitude','longitude','bands'], out_ds)}
        
    # Build Output xarray Dataset and assign data_vars array attributes
    out_xr = xr.Dataset(data_vars=data_vars, coords=coords, attrs= ds.attrs)
    # Assign Attributes from Original Datasets
    out_xr[var].attrs = ds[var].attrs
    #out_xr.coords['spatial_ref'] = ds.spatial_ref
    out_xr.coords['latitude'].attrs = loc['lon'].attrs
    out_xr.coords['longitude'].attrs = loc['lat'].attrs
    out_xr.rio.write_crs(ds.spatial_ref,inplace=True)
    #out_elev.attrs = loc['elev'].attrs
    # Mask Fill Values
    out_xr = out_xr.where(out_xr[var] != fill_value)
    
    return out_xr #, out_elev

def build_mask(filepath, quality_bands):
    """
    This function builds a single mask to apply based on the bands selected from an EMIT L2A Mask file.

    Parameters:
    filepath: an EMIT L2A Mask netCDF file.
    quality_bands: a list of bands (quality flags only) from the mask file that should be used in creation of  mask.

    Returns: 
    quality_mask: a numpy array that can be used with the emit_xarray function to apply a quality mask.
    """
    # Open Dataset
    mask_ds = xr.open_dataset(filepath,engine = 'h5netcdf')
    # Open Sensor band Group
    mask_parameters_ds = xr.open_dataset(filepath,engine = 'h5netcdf', group='sensor_band_parameters')
    # Print Flags used
    flags_used = mask_parameters_ds['mask_bands'].data[quality_bands]
    print(f'Flags used: {flags_used}')
    # Check for data bands and build mask
    if any(x in quality_bands for x in [5,6]):
        err_str = f'Selected flags include a data band (5 or 6) not just flag bands'
        raise AttributeError(err_str)
    else:
        quality_mask = np.sum(mask_ds['mask'][:,:,quality_bands].values, axis=-1)
        quality_mask[quality_mask > 1] = 1
    return(quality_mask)

#######################################################################################
#### Below here is under development/unused 
#######################################################################################
def point_extract(point_df,xarray_dataset):
    """
    This function extracts the lat and lon coordinates of a pandas dataframe from an EMIT xarray.Dataset using a 'nearest' method.

    Parameters:
    point_df: a pandas dataframe containing 'Latitude' and 'Longitude' columns.
    xarray_dataset: and EMIT dataset read in as an xarray.Dataset using the emit_xarray function.

    Returns: 
    out_df: a pandas dataframe containing columns from the input dataframe and the extracted variable values (i.e. reflectance).
    """
    # This could be improved by allowing for different variations of lat/lon

    # Select Pixels closest to Lat/Lon coordinates Provided and convert to dataframe.
    ex_df = xarray_dataset.sel(latitude=point_df.Latitude.to_xarray(),longitude=point_df.Longitude.to_xarray(),
                                         method='nearest').to_dataframe()
    # Pivot to wide-form of variable and create column headings from wavelengths
    ex_df = ex_df.pivot_tab+le(index=['index'],columns=['wavelengths'],values=['reflectance'])
    # Rename/Flatten Column Names
    colnames = []
    for col in ex_df.columns.values.tolist():
        name = '_'.join(map(str,col))
        colnames.append(name)
    ex_df.columns = colnames
    # Join original dataframe with extracted values
    out_df = point_df.join(ex_df)
    return(out_df)

def area_extract(area_gpd,xarray_dataset):
    """
    This function extracts the area of a single polygon in a geopandas dataframe from an EMIT xarray.Dataset.

    Parameters:
    point_df: a pandas dataframe containing 'Latitude' and 'Longitude' columns.
    xarray_dataset: and EMIT dataset read in as an xarray.Dataset using the emit_xarray function.

    Returns: 
    out_df: a pandas dataframe containing columns from the input dataframe and the extracted variable values (i.e. reflectance).
    """