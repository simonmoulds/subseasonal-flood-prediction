# Main entrypoint of the workflow.
# Please follow the best practices:
# https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html,
# in particular regarding the standardized folder structure mentioned there.

import platform
import math
import numpy as np
import pandas as pd
import pyarrow

from datetime import datetime, timedelta
from itertools import product

current_platform = platform.system()

# Discharge stations to include in prediction
METADATA = pd.read_parquet('resources/stations_selected.parquet')
CAMELS_STATIONS = list(
    METADATA['id'][(METADATA['has_camels'])].astype(int).astype(str)
)
STATIONS = CAMELS_STATIONS

c3s_start_year = 1981
c3s_end_year = 2016

C3S_VARIABLE = ['tp', 't2m']
C3S_YEARS = [yr for yr in range(c3s_start_year, c3s_end_year + 1)]
C3S_INIT_MONTH = [m for m in range(1, 13)]
C3S_INIT_TIME = [str(year) + str(month).zfill(2) + '01' for year, month in product(C3S_YEARS, C3S_INIT_MONTH)]

# Iterate over each year in the range
def haduk_time_ranges(start_year, end_year):
    rngs = []
    for year in range(start_year, end_year + 1):
        # Iterate over each month in the year
        for month in range(1, 13):
            # Calculate the first day of the month
            start_date = datetime(year, month, 1)
            # Calculate the last day of the month
            if month == 12:
                end_date = datetime(year + 1, 1, 1) - timedelta(days=1)
            else:
                end_date = datetime(year, month + 1, 1) - timedelta(days=1)

            start_date = start_date.strftime('%Y%m%d').strip()
            end_date = end_date.strftime('%Y%m%d').strip()

            rngs.append(start_date + '-' + end_date)

    return rngs

haduk_start_year = 1960
haduk_end_year = 2016

HADUK_VARIABLE = ['rainfall', 'tasmin', 'tasmax']
HADUK_YEARS = [yr for yr in range(haduk_start_year, haduk_end_year + 1)]
HADUK_TIME_RANGE = haduk_time_ranges(haduk_start_year, haduk_end_year)

# Modelling
# TODO put in config
PREDICTAND = ['Qd_mean']

# EFAS inputs
EFAS_YEARS = [yr for yr in range(1991, 2023)]
EFAS_MONTHS = [m for m in range(1, 13)]
EFAS_TIME = [str(year) + str(month).zfill(2) + '01' for year, month in product(EFAS_YEARS, EFAS_MONTHS)]

# Adjust if on dev platform
if current_platform == "Darwin":
    STATIONS = STATIONS[slice(0, 10)]
    C3S_VARIABLE = ['tp']
    C3S_YEARS = [yr for yr in range(2010, 2011)]
    C3S_INIT_MONTH = [m for m in range(1, 13)]
    C3S_INIT_TIME = [str(year) + str(month).zfill(2) + '01' for year, month in product(C3S_YEARS, C3S_INIT_MONTH)]

    HADUK_VARIABLE = ['rainfall']
    HADUK_YEARS = [yr for yr in range(2022, 2023)]
    HADUK_TIME_RANGE = haduk_time_ranges(2022, 2022)

    EFAS_YEARS = [1991, 2020]
    EFAS_MONTHS = [m for m in range(1, 13)]
    EFAS_TIME = [str(year) + str(month).zfill(2) + '01' for year, month in product(EFAS_YEARS, EFAS_MONTHS)]

EFAS_METADATA = pd.read_parquet('resources/efas_sites.parquet')
EFAS_METADATA = EFAS_METADATA[EFAS_METADATA['nrfa_id'].isin(STATIONS)]
EFAS_STATIONS = list(EFAS_METADATA['nrfa_id'].astype(str))

# Apply wildcard constraints
wildcard_constraints:
    station='|'.join([re.escape(x) for x in STATIONS]),
    efas_station='|'.join([re.escape(str(x)) for x in EFAS_STATIONS]),
    efas_time='|'.join([re.escape(str(x)) for x in EFAS_TIME]),
    c3s_variable='|'.join([re.escape(x) for x in C3S_VARIABLE]),
    c3s_init_time='|'.join([re.escape(str(x)) for x in C3S_INIT_TIME]),
    haduk_time_range='|'.join([re.escape(str(x)) for x in HADUK_TIME_RANGE]),
    predictand='|'.join([re.escape(x) for x in PREDICTAND])
