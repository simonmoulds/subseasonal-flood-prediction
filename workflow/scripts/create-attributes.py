#!/usr/bin/env python

import os
import numpy as np
import pandas as pd 

from pathlib import Path
from functools import reduce

stations = snakemake.params[0]
attr_output_file = Path(snakemake.output[0])
static_attributes = snakemake.config['prediction']['static_attributes']

# # TESTING:
# import yaml
# # read the yaml file
# with open('config/config.yml', 'r') as file:
#     cfg = yaml.safe_load(file)
# static_attributes = cfg['prediction']['static_attributes']

def clean_column(value):
    if isinstance(value, float) and value.is_integer():
        return str(int(value))  # Convert to int first, then to string
    return str(value)

# Static catchment attributes
attributes = pd.read_parquet("resources/stations_selected.parquet")
attributes['id'] = attributes['id'].apply(clean_column)
# attributes['id'] = attributes['id'].astype(str)
# attributes['id'] = attributes['id'].astype(int).astype(str)
attributes.set_index('id', inplace=True)
# attributes = attributes.loc[stations]
attributes = attributes.reset_index().rename({'id': 'gauge_id'}, axis=1)

# Filter stations with NaN values 
attributes = attributes.dropna(subset=static_attributes)

exclude = ['23001', '25001', '25009', '27090', '33051', '33066', '34002', '34012', '39130', '40004', '41005', '46008', '54011', '54038', '6011', '69020', '73013', '74008', '75005', '76011']

attributes = attributes[~attributes['gauge_id'].isin(exclude)]

attributes.to_csv(attr_output_file, index=False)
