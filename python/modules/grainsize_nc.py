# Predict grain size for an EMIT scene using the UCLA machine learning model
# Written Philip G. Brodrick
# Archived by David R. Thompson

# Modified from https://github.com/emit-sds/emit-sds-l3/blob/488ebd627cca6395a9d1d9ffe6cf3b6b42afc223/grainsize.py
# works with .nc rather than binary/hdr 
# Erik Bolch, ebolch@contractor.usgs.gov

# Required Model: https://zenodo.org/records/16700296

# Environment
# mamba env create -n grainsize_nc -c conda-forge python=3.10 scikit-learn=1.1.3 numpy=1.23.4 scipy=1.9.3 netcdf4 gdal spectral ipykernel

# Example
# python python\modules\grainsize_nc.py data\mineral_tutorial\EMIT_L2A_RFL_001_20230622T193324_2317313_010.nc data\multiRF.sav data\mineral_tutorial\EMIT_GRS_20230622T193324.nc

import argparse
import numpy as np
import netCDF4 as nc
from scipy import interpolate
import logging
import pickle
import os

def write_output(outdat, outfile, shape, band_names, rfl_file):
    logging.info(f'Writing output file {outfile}')
    # write in BIL
    with nc.Dataset(outfile, 'w', format='NETCDF4') as ds_out:
        ds_out.description = 'Predicted grainsize from RandomForestRegressor model'
        ds_out.createDimension('band', shape[2])
        ds_out.createDimension('downtrack', shape[0])
        ds_out.createDimension('crosstrack', shape[1])
        
        bands = ds_out.createVariable('band', str, ('band',))
        bands[:] = np.array(band_names)
        
        predictions = ds_out.createVariable('grainsize_predictions', np.float32, ('band', 'downtrack', 'crosstrack'), fill_value=-9999)
        predictions.units = 'unitless'
        predictions[:] = np.transpose(outdat, (2,0,1))

        # Add Location Group
        out_loc = ds_out.createGroup('location')

        with nc.Dataset(rfl_file, 'r') as src:
            loc_grp = src.groups['location']

            # Copy Dims
            needed_dims = set()
            for v in loc_grp.variables.values():
                needed_dims.update(v.dimensions)
            for dn in needed_dims:
                if dn in ds_out.dimensions:
                    continue
                if dn in loc_grp.dimensions:
                    d = loc_grp.dimensions[dn]
                elif dn in src.dimensions:
                    d = src.dimensions[dn]
                else:
                    raise ValueError(f"Dimension {dn} not found in source file")
                ds_out.createDimension(dn, None if d.isunlimited() else len(d))
            # Copy Group Attrs
            for name in loc_grp.ncattrs():
                out_loc.setattr(name, loc_grp.getncattr(name))
            
            # Copy Vars
            for var_name, var in loc_grp.variables.items():
                out_var = out_loc.createVariable(
                    var_name, var.datatype, var.dimensions
                )
                out_var.setncatts({k: var.getncattr(k) for k in var.ncattrs()})
                out_var[:] = var[:]

def spectral_derivative(rfl):
    # Everything in between
    d = [np.array(-1.5*rfl[:,:,0]+2*rfl[:,:,1]-0.5*rfl[:,:,2])]
    d.append(np.empty(d[0].shape))
    for i in range(1, 283):
      d.append((rfl[:,:,i+1]-rfl[:,:,i-1])/2)
    d.append(np.array(1.5*rfl[:,:,-1]-2*rfl[:,:,-2]+0.5*rfl[:,:,-3]))

    d = np.stack(d,axis=-1)
    return d

def main():

    parser = argparse.ArgumentParser(description="Translate to Rrs. and/or apply masks")
    parser.add_argument('rfl_file', type=str, metavar='OUTPUT')
    parser.add_argument('model_file', type=str, metavar='Band Depth file.  4 bands (G1 BD, G1 Ref, G2 BD, G2 Ref)')
    parser.add_argument('output_file', type=str)
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--log_level', type=str, default='INFO')
    args = parser.parse_args()

    if os.path.isfile(args.output_file):
        print('already found output file...terminating')
        exit()

    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    logging.info(f'Opening {args.rfl_file}')
    with nc.Dataset(args.rfl_file, 'r') as ds:
        logging.info(f'Loading Data...')
        rfl = ds['reflectance'][:]
        rfl_shape = rfl.shape
        wl = ds['sensor_band_parameters']['wavelengths'][:]

    logging.info(f'...Data Loaded')
    logging.info(f'Calculate Spectral Derivative')
    rfl_d = spectral_derivative(rfl)

    valid_wl = np.array([x > 454 and not (x > 1298 and x < 1505) and not (x > 1775 and x < 1980) and x < 2301 for x in wl])

    rfl = rfl[...,valid_wl]
    rfl_d = rfl_d[...,valid_wl]

    # rfl = np.append(rfl_d,rfl,axis=2)
    # rfl = rfl.reshape((rfl.shape[0]*rfl.shape[1],rfl.shape[2]))
    # rfl = np.nan_to_num(rfl)
    nrow, ncol, nb = rfl.shape
    rfl_merged = np.empty((nrow*ncol, nb*2), dtype=np.float32)
    rfl_merged[:,:nb] = rfl_d.reshape(-1, nb)
    rfl_merged[:,nb:] = rfl.reshape(-1,nb)
    rfl_merged = np.nan_to_num(rfl_merged, copy=False)

    del rfl
    del rfl_d

    model = pickle.load(open(args.model_file, 'rb'))
    pred = model.predict(rfl_merged)

    # normalize
    pred = pred / np.sum(pred,axis=-1)[:,np.newaxis]
    predc = np.cumsum(pred,axis=-1)

    size_classes = np.array([1500, 750, 375, 187.5, 93.75, 30, 1])
    ifuns = [interpolate.interp1d(predc[n,:],size_classes) if rfl_merged[n,-100] > 0 else np.array([-1]) for n in range(pred.shape[0])]
    del predc
    median_size = np.array([fun([0.5]) if fun != -1 else np.array([-1]) for fun in ifuns])

    pred = np.hstack([pred,median_size])

    # mask
    pred[rfl_merged[...,-100] <= 0,:] = -9999

    pred = pred.reshape((rfl_shape[0],rfl_shape[1], pred.shape[-1]))
    write_output(pred, args.output_file, pred.shape, ['S1','S2','S3','S4','S5','TSI','Clay','Median Grainsize'], args.rfl_file)


if __name__ == "__main__":
    main()