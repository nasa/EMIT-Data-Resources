# Run mineral abundance model
# Written Philip G. Brodrick
# Archived by David R. Thompson

# Modified from https://github.com/emit-sds/emit-sds-l3/blob/488ebd627cca6395a9d1d9ffe6cf3b6b42afc223/m4p.jl
# works with .nc files rather than binary/hdr 
# Erik Bolch, ebolch@contractor.usgs.gov

import numpy as np
import argparse
import logging
import xarray as xr
import pandas as pd
from functools import partial
import dask.distributed

def filter_by_names(in_df, mineral_names):
    out_mat = in_df[mineral_names].fillna(0).to_numpy()
    out_mat[out_mat==-1] = 1
    return out_mat

def mean_optical_path_length(absorption_coefficient, band_depth):
    BD = 1.0 - band_depth
    delta_k = absorption_coefficient
    MOPL = np.maximum(np.log(BD) / (-1.0 * delta_k), 1e-12)
    return MOPL

def model_3(mineral_groupings, delta_k, grain_diameter, density, band_depth):

    abs_abundance = np.zeros(len(delta_k), dtype=np.float64)
    m_return = np.zeros(len(delta_k), dtype=np.float64)
    T_return = np.zeros(len(delta_k), dtype=np.float64)

    # Group 1 Minerals
    if band_depth[1] > 0:
        # Adjust for Indexing (Mineral Index is 1 based)
        group_idx = int(band_depth[1]-1)
        present_in = np.nonzero(mineral_groupings[group_idx,:])[0]
        if len(present_in) > 0:
            MOPL = mean_optical_path_length(delta_k[present_in],
                                            band_depth[0])
            D = grain_diameter[present_in]/10000
            m = MOPL / (2*D)
            T = np.pi / 4 * D * np.sqrt(m/2)
            abs_abundance[present_in] = 3.0 * T * 100.0**2 * density[present_in]
            m_return[present_in] = m
            T_return[present_in] = T
    
    # Group 2 Minerals
    if band_depth[3] > 0:
        # Adjust for Indexing (Mineral Index is 1 based)
        group_idx = int(band_depth[3] -1)
        present_in = np.nonzero(mineral_groupings[group_idx,:])[0]
        if len (present_in) > 0:
            MOPL = mean_optical_path_length(delta_k[present_in],
                                            band_depth[2])
            D = grain_diameter[present_in]/10000
            m = MOPL / (2*D)
            T = np.pi / 4 * D * np.sqrt(m/2)
            abs_abundance[present_in] = 3.0 * T * 100.0**2 * density[present_in]
            m_return[present_in] = m
            T_return[present_in] = T

    mean_m = np.mean(m_return[m_return != 0]) if np.any(m_return != 0 ) else 0.0

    return abs_abundance, mean_m, m_return, T_return

def process_pixel(band_depth, quartz_size, quartz_massfrac, mineral_groupings, delta_k, grain_diameter, density, 
                    quartz_size_microns, average_mingrains, quartz_gs_lb, q_dens, model_style):
    """
    Process a single pixel (designed for dask/apply_ufunc)
    """
    m3, m3m, m_return, T_return = model_3(mineral_groupings, delta_k, grain_diameter, density, band_depth)

    if quartz_size is None:
        qs = quartz_size_microns
    else:
        qs = np.maximum(quartz_size, quartz_gs_lb)
    
    out_data = np.zeros(len(delta_k) + 1, dtype=np.float32)
    out_m = np.zeros(len(delta_k) + 1, dtype=np.float32)
    out_T = np.zeros(len(delta_k) + 1, dtype=np.float32)

    if model_style == "m3":
        out_data[:-1] = m3
        if quartz_massfrac is not None:
            non_quartz_mass = np.sum(m3)
            quartz_mass = non_quartz_mass * quartz_massfrac / (1 - quartz_massfrac)
            out_data[-1] = quartz_mass

        else:
            T_q = np.pi / 4 * (qs / 10000.0) / average_mingrains * np.sqrt(m3m / 2)
            qmass = 3 * T_q * 100**2 * q_dens
            out_data[-1] = qmass
    else:
        raise ValueError(f"Unsuppoerted model_style: {model_style}")

    out_m[:-1] = m_return.astype(np.float32)
    out_m[-1] = np.float32(m3m)
    out_T[:-1] = T_return.astype(np.float32)
    return out_data, out_m, out_T

def main():
    parser = argparse.ArgumentParser(description="Estimate mineral abundance.")
    parser.add_argument('output_base', type=str, metavar='OUTPUT')
    parser.add_argument('band_depth_file', type=str, metavar='Band Depth file.  4 bands (G1 BD, G1 Ref, G2 BD, G2 Ref)')
    parser.add_argument('quartz_size_microns', type=np.float64)
    parser.add_argument('average_mingrains', type=np.float64)
    parser.add_argument('--quartz_size_file', type=str, default=None)
    parser.add_argument('--quartz_file_scaling_factor', type=np.float64, default=1)
    parser.add_argument('--quartz_uncertainty_delta', type=np.float64, default=0)
    parser.add_argument('--quartz_gs_lb', type=np.float64, default=0)
    parser.add_argument('--quartz_massfrac_file', type=str, default=None)
    parser.add_argument('--band_depth_unc_file', type=str, default=None)
    parser.add_argument('--model_style', type=str, default="m3")
    parser.add_argument('--q_dens', type=np.float64, default=2.65)
    parser.add_argument('--mineral_groupings_matrix', type=str, default='data/mineral_grouping_matrix_20231113.csv')
    parser.add_argument('--abundance_metadata', type=str, default='data/abundance_metadata20240123.csv')
    parser.add_argument('--sum_grainsize', action="store_true")
    parser.add_argument('--output_m', action="store_true")
    parser.add_argument('--log_file', type=str, default=None)
    parser.add_argument('--log_level', type=str, default='INFO')
    # The dask defaults work well for openscapes 2i2c
    parser.add_argument('--n_workers', type=int, default=3, help='Number of Dask workers (defaults to CPU count if processes=True)')
    parser.add_argument('--threads_per_worker', type=int, default=1, help='Threads per Dask worker')
    parser.add_argument('--processes', action='store_true', help='Use processes instead of threads for Dask client')
    parser.add_argument('--chunk_size', type=float, default=256, help='Target crosstrack and downtrack chunk size')    
    args = parser.parse_args()
    
    if args.log_file is None:
        logging.basicConfig(format='%(message)s', level=args.log_level)
    else:
        logging.basicConfig(format='%(message)s', level=args.log_level, filename=args.log_file)

    # Initialize Dask Client
    client = dask.distributed.Client(
        n_workers=args.n_workers,
        threads_per_worker=args.threads_per_worker,
        processes=args.processes,
        memory_limit='2GB'
    )

    logging.info(f"Dask client started: {client}")

    # Read band depth data
    band_depth_dt = xr.open_datatree(args.band_depth_file)
    band_depth = band_depth_dt['/'].to_dataset().to_array("band").transpose(
        "downtrack", "crosstrack","band").astype(np.float64).chunk(
            {"downtrack":args.chunk_size,"crosstrack":args.chunk_size,"band":-1})
    band_depth = band_depth.persist()
    
    #TODO Read band depth uncertainty
    band_depth_unc = None
    if args.band_depth_unc_file is not None:
        logging.warning("Band depth uncertainty file provided but not yet implemented. Ignoring.")

    #TODO Read in Quartz data
    if args.quartz_massfrac_file and args.quartz_size_file:
        logging.error("Can only have one of quartz_massfrac_file or quartz_size_file")
        exit(1)

    quartz_massfrac = None
    if args.quartz_massfrac_file:
        logging.warning("Quartz mass frac file provided but not yet implemented. Ignoring.")

    quartz_size = None
    if args.quartz_size_file:
       quartz_size_ds = xr.open_dataset(args.quartz_size_file)
       quartz_size = (quartz_size_ds['grainsize_predictions'].sel(band='Median Grainsize').astype(np.float64).chunk(
                {"downtrack":args.chunk_size,"crosstrack":args.chunk_size}) + args.quartz_uncertainty_delta) * args.quartz_file_scaling_factor
       quartz_size = quartz_size.persist()
       
    # Load metadata
    abundance_metadata = pd.read_csv(args.abundance_metadata)
    mineral_groupings = pd.read_csv(args.mineral_groupings_matrix)

    # Convert metadata to numpy arrays
    delta_k = abundance_metadata['delta_k_absorption_coefficient'].to_numpy().astype(np.float64)
    grain_diameter = abundance_metadata['grain_diameter_microns'].to_numpy().astype(np.float64)
    density = abundance_metadata['density_g_cc'].to_numpy().astype(np.float64)
    
    # Filter mineral groupings
    mi_header = mineral_groupings.columns
    mineral_names = [
        col for col in mi_header
        if mi_header.tolist().index('Calcite') <= mi_header.tolist().index(col) <= mi_header.tolist().index('Vermiculite')
        ]

    mineral_groupings = filter_by_names(mineral_groupings, mineral_names)
    num_minerals = len(mineral_names) + 1
    out_mineral_names = mineral_names + ["Quartz+Feldspar"]

    # Set up inputs for process_pixel function
    pixel_func = partial(process_pixel,
              quartz_massfrac=quartz_massfrac,
              mineral_groupings=mineral_groupings,
              delta_k=delta_k,
              grain_diameter=grain_diameter,
              density=density,
              quartz_size_microns=args.quartz_size_microns,
              average_mingrains=args.average_mingrains,
              quartz_gs_lb=args.quartz_gs_lb,
              q_dens=args.q_dens,
              model_style=args.model_style
              )
    
    scatter, light = None, None

    # Run Process Pixel Func
    abs_abundance, scatter, light = xr.apply_ufunc(
        pixel_func,
        band_depth,
        quartz_size,
        input_core_dims=[['band'], []],
        output_core_dims=[['minerals'], ['minerals'],['minerals']],
        vectorize=True,
        dask='parallelized',
        output_dtypes=[np.float32, np.float32, np.float32],
        dask_gufunc_kwargs={'allow_rechunk':True, 'output_sizes':{'minerals':num_minerals}}
        )
    
    # Add Fill Values and Calc Relative Abundance
    total_abundance = abs_abundance.sum(dim='minerals')
    abs_abundance= xr.where(total_abundance <= 0, -9999.0, abs_abundance)
    rel_abundance = abs_abundance / total_abundance.where(total_abundance > 0)       
        
    # Create Root Dataset
    ds = xr.Dataset(data_vars={"abs_abundance":abs_abundance,
                           "rel_abundance": rel_abundance})
    ds = ds.assign_attrs(band_depth.attrs)
    ds.attrs['title'] = "Locally computed EMIT L2B Estimated Mineral Abundance 60 m"

    if args.output_m:
        ds['scatter'] = scatter
        ds['light'] = light

    ds = ds.compute()

    # Set Root Group Encoding
    encoding = {
        "dtype":"float32",
        "_FillValue": np.float32(-9999),
        "zlib": True,
        "complevel":2,
        "shuffle":True
    }
    for var in ds.data_vars:
        ds[var].encoding.update(encoding)

    # Create Output DataTree
    dt = xr.DataTree(dataset=ds)

    # Copy Location from Band Depth DataTree
    dt['location'] = band_depth_dt['location'].copy()
    # Create Mineral Metadata Group and Add to DataTree
    dt['mineral_metadata'] = xr.Dataset({"name":(("minerals"),out_mineral_names)})

    # Write Output - Openable with emit_xarray
    dt.to_netcdf(args.output_base, engine='netcdf4', mode='w')
    client.close()

if __name__ == "__main__":
    main()






