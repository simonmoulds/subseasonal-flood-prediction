#!/usr/bin/env python3

import xarray
import numpy as np

# Parse Snakemake wildcards/parameters that are common to all rules
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# TESTING
# input_file = 'data/C3S_hindcast_daily_SEAS5/ecmf_SEAS5-v20171101_hindcast_S20120201_tp_daily.nc'
# outfile = 'data/test.nc'

x = xarray.open_dataset(input_file, mask_and_scale=True)
newvar = x['tp'].isel(time=slice(0, 90)) # Select first 90 days to save space
newvar = newvar.diff(dim='time')
newvar = xarray.concat([x['tp'].isel(time=[0]), newvar], dim='time')
newvar = np.maximum(newvar, 0.)
x['tp'] = newvar
x.to_netcdf(output_file)
