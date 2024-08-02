#!/usr/bin/env python

import os
import numpy as np
import pandas as pd 
import xarray as xr

from pathlib import Path
from functools import reduce

# This script creates pseudo-forecasts for training the model

station_id = str(snakemake.wildcards['station'])
member = str(snakemake.wildcards['member'])
streamflow_input_file = Path(snakemake.input[0])
hindcast_input_file = Path(snakemake.input[1])
forecast_input_file = Path(snakemake.input[2])
seq_length = int(snakemake.config['prediction']['seq_length'])
forecast_seq_length = int(snakemake.config['prediction']['forecast_seq_length'])
ts_output_file = Path(snakemake.output[0])

# # Testing 
# station_id = '97002'
# member = '0'
# streamflow_input_file = 'results/preprocessing/streamflow/timeseries/daily/97002.parquet'
# hindcast_input_file = 'results/preprocessing/HadUK/97002.parquet'
# forecast_input_file = 'results/preprocessing/C3S/97002_bc.parquet'
# seq_length = int(365)
# forecast_seq_length = int(45)

# Dynamic input data
q_input = pd.read_parquet(streamflow_input_file)
q_input['time'] = pd.to_datetime(q_input['time'])
q_input = q_input.rename({'time': 'date'}, axis=1)

hcst_input = pd.read_parquet(hindcast_input_file)
hcst_input['time'] = pd.to_datetime(hcst_input['time'])
hcst_input = hcst_input.pivot(columns='variable', index=['id', 'time'], values='value')
hcst_input = hcst_input.reset_index()
hcst_input = hcst_input.rename({'time': 'date'}, axis=1)
hcst_input['tmean'] = hcst_input[['tasmin', 'tasmax']].mean(axis=1)

fcst_input = pd.read_parquet(forecast_input_file)
fcst_input['member'] = fcst_input['member'].astype(int).astype(str)
fcst_input['init_time'] = pd.to_datetime(fcst_input['init_time'])
fcst_input['time'] = pd.to_datetime(fcst_input['time'])
fcst_input['lead_time'] = fcst_input['lead_time'].astype(int)
fcst_input = fcst_input[fcst_input['member'] == member]
index_vars = ['id', 'system', 'member', 'init_time', 'lead_time', 'time', 'lead_time_month', 'init_year']
fcst_input = fcst_input.pivot(columns='variable', index=index_vars, values='value')
fcst_input = fcst_input.reset_index()
fcst_input = fcst_input.rename({'t2m': 'tmean_forecast', 'tp': 'rainfall_forecast'}, axis=1)
fcst_input['tmean_forecast'] = fcst_input['tmean_forecast'] - 273.15 

# Create an input sequence for each initialization time
init_times = pd.to_datetime(fcst_input['init_time'].unique())

# # TESTING
# init_times = init_times[0:10]

hcst_input_list = []
fcst_input_list = []
for init_tm in init_times:
    # Hindcast component
    hcst_start_date = init_tm - pd.Timedelta(days=seq_length)
    hcst_end_date = init_tm
    hcst_ts = pd.DataFrame.from_dict({'date': pd.date_range(start=hcst_start_date, end=hcst_end_date, freq='1D')})
    hcst_input0 = hcst_ts.merge(hcst_input, on=['date'])
    hcst_input0 = hcst_input0.merge(q_input, on=['id', 'date'])
    hcst_input0.drop('id', axis=1, inplace=True)
    hcst_input0.set_index('date', inplace=True)
    hcst_input_list.append(hcst_input0)

    # Forecast component
    fcst_start_date = init_tm + pd.Timedelta(days=1)
    fcst_end_date = init_tm + pd.Timedelta(days=forecast_seq_length)
    fcst_ts = pd.DataFrame.from_dict({'time': pd.date_range(start=fcst_start_date, end=fcst_end_date, freq='1D')})
    fcst_input0 = fcst_input[fcst_input['init_time'] == init_tm]
    fcst_input0 = fcst_ts.merge(fcst_input0, on=['time'])
    
    # Remove superfluous columns, add index 
    fcst_input0 = fcst_input0.rename({'init_time': 'date'}, axis=1)
    fcst_input0.set_index(['date', 'lead_time'], inplace=True)
    fcst_input0 = fcst_input0[[col for col in fcst_input0.columns if col.endswith('_forecast')]]
    fcst_input_list.append(fcst_input0)

hcst_input = pd.concat(hcst_input_list)
hcst_input = hcst_input[~hcst_input.index.duplicated(keep='first')]
fcst_input = pd.concat(fcst_input_list)

ds_hcst = xr.Dataset.from_dataframe(hcst_input)
ds_fcst = xr.Dataset.from_dataframe(fcst_input)

ds = ds_fcst.merge(ds_hcst) 
ds.to_netcdf(ts_output_file)

# # Merge discharge and meteo data
# dyn_input = q_input.merge(met_input, on=['id', 'time'])

# # Add forecast data column 
# dyn_input['rainfall_forecast'] = dyn_input['rainfall']
# dyn_input['tmean_forecast'] = dyn_input['tmean']
# dyn_input = dyn_input.rename({'time': 'date'}, axis=1)
# dyn_input['date'] = pd.to_datetime(dyn_input['date'])

# # Ensure complete time series
# ts = pd.DataFrame.from_dict({'date': pd.date_range(start=train_start_date, end=train_end_date, freq='1D')})
# dyn_input = dyn_input.merge(ts, on=['date'], how='outer')
# dyn_input.set_index('date', inplace=True)

# dyn_input_fcst_list = []
# dyn_input_hcst_list = []
# for init_tm in hindcast_init_time: 
#     row_idx = dyn_input.index.get_loc(init_tm)

#     # Hindcast portion:
#     dyn_input_hcst = dyn_input.iloc[[row_idx]].reset_index()
#     dyn_input_hcst.set_index('date', inplace=True)
#     dyn_input_hcst = dyn_input_hcst[[col for col in dyn_input_hcst.columns if not col.endswith('_forecast')]]
#     dyn_input_hcst_list.append(dyn_input_hcst)

#     # Forecast portion
#     dyn_input_fcst = dyn_input.iloc[slice(row_idx, row_idx + hindcast_seq_length)].reset_index()
#     dyn_input_fcst['date'] = dyn_input_fcst['date'].values[0]
#     dyn_input_fcst['lead_time'] = [i for i in range(hindcast_seq_length)]
#     dyn_input_fcst.set_index(['date', 'lead_time'], inplace=True)
#     dyn_input_fcst = dyn_input_fcst[[col for col in dyn_input_fcst.columns if col.endswith('_forecast')]]
#     dyn_input_fcst_list.append(dyn_input_fcst)

# dyn_input_hcst = pd.concat(dyn_input_hcst_list)
# dyn_input_fcst = pd.concat(dyn_input_fcst_list)

# ds_hcst = xr.Dataset.from_dataframe(dyn_input_hcst)
# ds_fcst = xr.Dataset.from_dataframe(dyn_input_fcst)

# ds = ds_fcst.merge(ds_hcst) 
# ds.to_netcdf(ts_output_file)
# attributes.to_csv(attr_output_file)
