#!/usr/bin/env python

import pandas as pd
import numpy as np

from pathlib import Path

attributes_filename = snakemake.input[0]
stations = snakemake.params[0]
basins_output_file = Path(snakemake.output[0])

attributes = pd.read_csv(attributes_filename, dtype={'gauge_id': str})
attributes = attributes[attributes['gauge_id'].isin(stations)]

selected_stations = list(attributes['gauge_id'])

with open(basins_output_file, 'w') as file:
    for station in selected_stations:
        file.write(station + '\n')