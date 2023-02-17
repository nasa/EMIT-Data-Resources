"""
This Module has the functions related to working with an EMIT dataset. This includes doing things
like opening and flattening the data to work in xarray, orthorectification, and extracting point and area samples.

Author: Erik Bolch, ebolch@contractor.usgs.gov 

Last Updated: 02/16/2022

TO DO: 
- Fix elevation orthorectification - all values end up as -9999
- Investigate reducing memory usage
- Test/Improve flexibility for applying the GLT to modified datasets
"""
# Packages used
import netCDF4 as nc
from osgeo import gdal
import numpy as np
import math
import pandas as pd
import xarray as xr
import rasterio as rio

def emit_xarray(filepath, ortho=True, qmask=None, unpacked_bmask=None): 
    """
        This function utilizes other functions in this module to streamline opening an EMIT dataset as an xarray.Dataset.
        
        Parameters:
        filepath: a filepath to an EMIT netCDF file
        ortho: True or False, whether to orthorectify the dataset or leave in crosstrack/downtrack coordinates.
        qmask: a numpy array output from the quality_mask function used to mask pixels based on quality flags selected in that function. Any non-orthorectified array with the proper crosstrack and downtrack dimensions can also be used.
        unpacked_bmask: a numpy array from  the band_mask function that can be used to mask band-specific pixels that have been interpolated.
                        
        Returns:
        out_xr: an xarray.Dataset constructed based on the parameters provided.

        """
    # Read in Data as Xarray Datasets
    ds = xr.open_dataset(filepath,engine = 'h5netcdf')
    loc = xr.open_dataset(filepath, engine = 'h5netcdf', group='location')
    wvl = xr.open_dataset(filepath, engine = 'h5netcdf', group='sensor_band_parameters') 
    
    # Building Flat Dataset from Components
    data_vars = {**ds.variables} 
    coords = {'downtrack':(['downtrack'], ds.downtrack.data),'crosstrack':(['crosstrack'],ds.crosstrack.data), **loc.variables, **wvl.variables}
    out_xr = xr.Dataset(data_vars=data_vars, coords = coords, attrs= ds.attrs)
    
    # Apply Quality and Band Masks
    if qmask is not None:
        out_xr[list(out_xr.data_vars)[0]].values[qmask == 1] = np.nan
    if unpacked_bmask is not None:
        out_xr[list(out_xr.data_vars)[0]].values[unpacked_bmask == 1] = np.nan               
    
    if ortho is True:
       out_xr = ortho_xr(out_xr)
          
    return out_xr


# Function to Calculate the Lat and Lon Vectors/Coordinate Grid
def coord_vects(ds):
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
    dim_x = ds.glt_x.shape[1]
    dim_y = ds.glt_x.shape[0]
    lon = np.zeros(dim_x)
    lat = np.zeros(dim_y)
    # Note: no rotation for EMIT Data
    for x in np.arange(dim_x):
        x_geo = GT[0] + x * GT[1]
        lon[x] = x_geo
    for y in np.arange(dim_y):
        y_geo = GT[3] + y * GT[5]
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
    out_ds = np.full((glt_array.shape[0], glt_array.shape[1], ds_array.shape[-1]), fill_value, dtype=np.float32)
    valid_glt = np.all(glt_array != GLT_NODATA_VALUE, axis=-1)
    
    # Adjust for One based Index
    glt_array[valid_glt] -= 1 
    out_ds[valid_glt, :] = ds_array[glt_array[valid_glt, 1], glt_array[valid_glt, 0], :]
    return out_ds

def ortho_xr(ds, GLT_NODATA_VALUE=0, fill_value = -9999):
    """
    This function applies the GLT to variables within an EMIT dataset that has been read into the format provided in by the emit_xarray function.

    Parameters:
    ds: an xarray dataset produced by emit_xarray
    GLT_NODATA_VALUE: no data value for the GLT tables, 0 by default
    fill_value: the fill value for EMIT datasets, -9999 by default

    Returns:
    ortho_ds: an orthocorrected xarray dataset.  
    
    """
    # Build glt_ds

    glt_ds = np.nan_to_num(np.stack([ds['glt_x'].data,ds['glt_y'].data],axis=-1),nan=GLT_NODATA_VALUE).astype(int)  
    

    # List Variables
    var_list = list(ds.data_vars)
    
    # Create empty dictionary for orthocorrected data vars
    data_vars = {}   

    # Extract Rawspace Dataset Variable Values (Typically Reflectance)    
    for var in var_list:
        raw_ds = ds[var].data

        # Apply GLT to dataset
        out_ds = apply_glt(raw_ds,glt_ds)

        del raw_ds
        #Update variables
        data_vars[var] = (['latitude','longitude','bands'], out_ds)
    
    # Calculate Lat and Lon Vectors
    lon, lat = coord_vects(ds) # Reorder this function to make sense in case of multiple variables

    # Apply GLT to elevation
    #elev_ds = apply_glt(ds['elev'].data[:,:,np.newaxis],glt_ds)
    
    # Delete glt_ds - no longer needed
    del glt_ds
    
    # Create Coordinate Dictionary
    coords = {'latitude':(['latitude'],lat), 'longitude':(['longitude'],lon), **ds.coords}# unpack to add appropriate coordinates 
    
    # Remove Unnecessary Coords
    for key in ['downtrack','crosstrack','lat','lon','glt_x','glt_y','elev']:
        del coords[key]
    
    # Add Orthocorrected Elevation
    #coords['elev'] = (['latitude','longitude'], np.squeeze(elev_ds))

    # Build Output xarray Dataset and assign data_vars array attributes
    out_xr = xr.Dataset(data_vars=data_vars, coords=coords, attrs=ds.attrs)
    
    del out_ds
    # Assign Attributes from Original Datasets
    out_xr[var].attrs = ds[var].attrs
    out_xr.coords['latitude'].attrs = ds['lon'].attrs
    out_xr.coords['longitude'].attrs = ds['lat'].attrs
    
    # Add Spatial Reference in recognizable format
    out_xr.rio.write_crs(ds.spatial_ref,inplace=True)

    # Mask Fill Values
    out_xr[var].data[out_xr[var].data == fill_value] = np.nan
    
    return out_xr  

def quality_mask(filepath, quality_bands):
    """
    This function builds a single layer mask to apply based on the bands selected from an EMIT L2A Mask file.

    Parameters:
    filepath: an EMIT L2A Mask netCDF file.
    quality_bands: a list of bands (quality flags only) from the mask file that should be used in creation of  mask.

    Returns: 
    qmask: a numpy array that can be used with the emit_xarray function to apply a quality mask.
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
        qmask = np.sum(mask_ds['mask'][:,:,quality_bands].values, axis=-1)
        qmask[qmask > 1] = 1
    return(qmask)

def band_mask(filepath):
    """
    This function unpacks the packed band mask to apply to the dataset. Can be used manually or as an input in the emit_xarray() function.

    Parameters:
    filepath: an EMIT L2A Mask netCDF file.
    packed_bands: the 'packed_bands' dataset from the EMIT L2A Mask file.

    Returns: 
    band_mask: a numpy array that can be used with the emit_xarray function to apply a band mask.
    """
    # Open Dataset
    mask_ds = xr.open_dataset(filepath,engine = 'h5netcdf')
    # Open band_mask and convert to uint8
    bmask = mask_ds.band_mask.data.astype('uint8')
    # Print Flags used
    unpacked_bmask = np.unpackbits(bmask,axis=-1)
    # Remove bands > 285
    unpacked_bmask = unpacked_bmask[:,:,0:285]
    # Check for data bands and build mask
    return(unpacked_bmask)
