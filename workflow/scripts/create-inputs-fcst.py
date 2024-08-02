#!/usr/bin/env python

import os
import numpy as np
import pandas as pd 
import xarray as xr

from pathlib import Path
from functools import reduce

input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Attributes
CAMELS_GB_ATTRIBUTE_FILES = [
    'CAMELS_GB_climatic_attributes.csv',
    'CAMELS_GB_humaninfluence_attributes.csv',
    'CAMELS_GB_hydrogeology_attributes.csv',
    'CAMELS_GB_hydrologic_attributes.csv',
    'CAMELS_GB_hydrometry_attributes.csv',
    'CAMELS_GB_landcover_attributes.csv',
    'CAMELS_GB_soil_attributes.csv',
    'CAMELS_GB_topographic_attributes.csv'
]

# Forecast initialization times
start_date = '1994-01-01'
end_date = '2014-12-01'
FORECAST_INIT_TIME = pd.date_range(start=start_date, end=end_date, freq='MS')
FORECAST_N_DAYS = 45

# Ensemble forecast members 
MEMBERS = [i for i in range(0, 25)]

DATA_DIR = Path('data/camels_gb/8344e4f3-d2ea-44f5-8afa-86d2987543a9/data')

OUTPUT_DIR = Path('data/camels_gb_forecast')

def main():

    # Create directories to store output
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR / 'time_series', exist_ok=True)
    os.makedirs(OUTPUT_DIR / 'attributes', exist_ok=True)

    # Subselect IDs
    x = pd.read_csv(DATA_DIR / 'CAMELS_GB_climatic_attributes.csv')
    ids = x['gauge_id'].values
    ids_selected = x.sample(8, random_state=42)['gauge_id'].values
    basins_filename = 'data/8_basins_train_set.txt'
    np.savetxt(basins_filename, ids_selected, fmt='%i')

    # Load attributes and convert to single dataframe
    attributes = []
    for fn in CAMELS_GB_ATTRIBUTE_FILES:
        attr = pd.read_csv(os.path.join(DATA_DIR, fn))
        attributes.append(attr)

    attributes = reduce(lambda x, y: pd.merge(x, y, on='gauge_id'), attributes)

    # Create pseudo-forecasts for experimentation
    test_sets = {}
    for basin in ids_selected:
        fname = f'CAMELS_GB_hydromet_timeseries_{basin}_19701001-20150930.csv'
        dyn_input = pd.read_csv(os.path.join(DATA_DIR, 'timeseries', fname), parse_dates=[0])

        # Add forecast data column 
        dyn_input['precipitation_forecast'] = dyn_input['precipitation']
        dyn_input['temperature_forecast'] = dyn_input['temperature']

        # Write training data 
        dyn_input.set_index('date', inplace=True)
        ds = xr.Dataset.from_dataframe(dyn_input)

        pseudo_basins = []
        for init_tm in FORECAST_INIT_TIME:
            hindcast_start_time = init_tm - pd.Timedelta(days=365)
            forecast_end_time = init_tm + pd.Timedelta(days=FORECAST_N_DAYS)
            dyn_input0 = dyn_input[(dyn_input.index >= hindcast_start_time) & (dyn_input.index <= forecast_end_time)].copy()

            forecast_cols = [col for col in dyn_input0.columns if col.endswith('forecast')]
            dyn_input0.loc[:init_tm, forecast_cols] = np.nan 

            ds0 = xr.Dataset.from_dataframe(dyn_input0)

            pseudo_basin_name = str(basin) + '_S' + init_tm.strftime('%d%m%Y')
            pseudo_basins.append(pseudo_basin_name)

            attr_basin = attributes[attributes['gauge_id'] == basin].copy()
            attr_basin['gauge_id'] = pseudo_basin_name
            attributes = pd.concat([attributes, attr_basin]).reset_index(drop=True)

            # Write output data 
            ds0.to_netcdf(os.path.join(OUTPUT_DIR, 'time_series', pseudo_basin_name + '.nc'))

        test_sets[basin] = pseudo_basins 

        ds.to_netcdf(os.path.join(OUTPUT_DIR, 'time_series', str(basin) + '.nc'))

    attributes.set_index('gauge_id', inplace=True)
    attributes.to_csv(os.path.join(OUTPUT_DIR, 'attributes', 'attributes.csv'))

    # Write test sets 
    testsets = pd.DataFrame.from_dict(test_sets)
    testsets['init_time'] = FORECAST_INIT_TIME 
    testsets.set_index('init_time', inplace=True)
    for init_tm in FORECAST_INIT_TIME:
        basins = testsets[testsets.index == init_tm]
        if basins.shape[0] == 1: 
            basins = basins.values.flatten()  
            basins_filename = 'data/8_basins_test_set_S' + init_tm.strftime('%d%m%Y') + '.txt'
            np.savetxt(basins_filename, basins, fmt='%s')

if __name__ == '__main__':
    main()

