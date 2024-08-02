#!/usr/bin/env python

import os
import numpy as np
import pandas as pd 
import xarray as xr

from pathlib import Path
from functools import reduce

# This script creates pseudo-forecasts for training the model

station_id = str(snakemake.wildcards['station'])
streamflow_input_file = Path(snakemake.input[0])
climate_input_file = Path(snakemake.input[1])
ts_output_file = Path(snakemake.output[0])
# attr_output_file = Path(snakemake.output[1])

train_start_date = pd.Timestamp(snakemake.config['prediction']['start_date'])
train_end_date = pd.Timestamp(snakemake.config['prediction']['end_date'])
hindcast_seq_length = int(snakemake.config['prediction']['forecast_seq_length'])

# # Testing 
# station_id = '97002'
# streamflow_input_file = 'results/preprocessing/streamflow/timeseries/daily/97002.parquet'
# climate_input_file = 'results/preprocessing/HadUK/97002.parquet'
# train_start_date = pd.Timestamp('1990/01/01')
# train_end_date = pd.Timestamp('2004/12/31')
# hindcast_seq_length = int(45)

hindcast_init_time = pd.date_range(start=train_start_date, end=train_end_date, freq='1D')
hindcast_init_time = hindcast_init_time[:-(hindcast_seq_length)]

# # Static catchment attributes
# attributes = pd.read_parquet("resources/stations_selected.parquet")
# attributes['id'] = attributes['id'].astype(int).astype(str)
# attributes.set_index('id', inplace=True)
# attributes = attributes.loc[[station_id]]
# attributes = attributes.reset_index().rename({'id': 'gauge_id'}, axis=1)

# Dynamic input data
q_input = pd.read_parquet(streamflow_input_file)
met_input = pd.read_parquet(climate_input_file)
met_input = met_input.pivot(columns='variable', index=['id', 'time'], values='value')
met_input = met_input.reset_index()
# NOTE HadUK does not give mean temp, so estimate it here
# NOTE In the future we will use ERA5 data 
met_input['tmean'] = met_input[['tasmax', 'tasmin']].mean(axis=1)

# Merge discharge and meteo data
dyn_input = q_input.merge(met_input, on=['id', 'time'])

# Add forecast data column 
dyn_input['rainfall_forecast'] = dyn_input['rainfall']
dyn_input['tmean_forecast'] = dyn_input['tmean']
dyn_input = dyn_input.rename({'time': 'date'}, axis=1)
dyn_input['date'] = pd.to_datetime(dyn_input['date'])

# Ensure complete time series
ts = pd.DataFrame.from_dict({'date': pd.date_range(start=train_start_date, end=train_end_date, freq='1D')})
dyn_input = dyn_input.merge(ts, on=['date'], how='outer')
dyn_input.set_index('date', inplace=True)

dyn_input_fcst_list = []
dyn_input_hcst_list = []
for init_tm in hindcast_init_time: 
    row_idx = dyn_input.index.get_loc(init_tm)

    # Hindcast portion:
    dyn_input_hcst = dyn_input.iloc[[row_idx]].reset_index()
    dyn_input_hcst.set_index('date', inplace=True)
    dyn_input_hcst = dyn_input_hcst[[col for col in dyn_input_hcst.columns if not col.endswith('_forecast')]]
    dyn_input_hcst.drop('id', axis=1, inplace=True)
    dyn_input_hcst_list.append(dyn_input_hcst)

    # Forecast portion 
    # +1 because seasonal forecasts start one day after the nominal issue date
    # (i.e. for a forecast issued 01/01/2010, the first forecast will be 02/01/2010)
    dyn_input_fcst = dyn_input.iloc[slice(row_idx + 1, row_idx + 1 + hindcast_seq_length)].reset_index()
    dyn_input_fcst['date'] = dyn_input_fcst['date'].values[0]
    dyn_input_fcst['lead_time'] = [(i + 1) for i in range(hindcast_seq_length)]
    dyn_input_fcst.set_index(['date', 'lead_time'], inplace=True)
    dyn_input_fcst = dyn_input_fcst[[col for col in dyn_input_fcst.columns if col.endswith('_forecast')]]
    dyn_input_fcst_list.append(dyn_input_fcst)

dyn_input_hcst = pd.concat(dyn_input_hcst_list)
dyn_input_fcst = pd.concat(dyn_input_fcst_list)

ds_hcst = xr.Dataset.from_dataframe(dyn_input_hcst)
ds_fcst = xr.Dataset.from_dataframe(dyn_input_fcst)

ds = ds_fcst.merge(ds_hcst) 
ds.to_netcdf(ts_output_file)
# attributes.to_csv(attr_output_file)
