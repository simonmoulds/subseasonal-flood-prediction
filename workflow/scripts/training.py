# import os 
# import pandas as pd
# import xarray 
# import pickle
# import matplotlib.pyplot as plt
# import importlib 

from pathlib import Path

# from neuralhydrology.training.basetrainer import BaseTrainer
from neuralhydrology.training.train import start_training
# from neuralhydrology.evaluation import get_tester
# from neuralhydrology.evaluation.evaluate import start_evaluation 
from neuralhydrology.utils.config import Config

# Snakemake stuff
basin_file = Path(snakemake.input[0])
base_config_file = Path(snakemake.params[0])
output_file = Path(snakemake.output[0])

# # Testing:
# basin_file = Path('./results/preprocessing/merged/basin_list_8.txt')
# base_config_file = Path('./resources/config.yml.example')
# output_file = Path('./results/training/latest_run_dir.txt')
run_dir = output_file.parent

cfg = Config(base_config_file)
cfg_dict = {
    'experiment_name': 'subseasonal_forecast_run',
    'run_dir': run_dir,
    'train_basin_file': basin_file,
    'validation_basin_file': basin_file,
    'test_basin_file': basin_file,
    # training, validation and test time periods (format = 'dd/mm/yyyy')
    # 'train_start_date': '01/10/1999',
    # 'train_end_date': '30/09/2008',
    # 'validation_start_date': '01/10/2008',
    # 'validation_end_date': '30/09/2014',
    # 'test_start_date': '01/10/2008',
    # 'test_end_date': '30/09/2014',
    'train_start_date': '01/10/1993',
    'train_end_date': '30/09/2003',
    'validation_start_date': '01/10/1993',
    'validation_end_date': '30/09/2003',
    'test_start_date': '01/10/2003',
    'test_end_date': '30/09/2016',
    # fixed seed, leave empty to use a random seed
    'seed': 42,
    # which GPU (id) to use [in format of cuda:0, cuda:1 etc, or cpu or None]
    'device': 'cuda:1',
    # base model type [cudalstm, customlstm, ealstm, embcudalstm, mtslstm, gru, transformer]
    # (has to match the if statement in modelzoo/__init__.py)
    'model': 'handoff_forecast_lstm',
    # Defines which time steps are used to calculate the loss. Can't be larger than seq_length.
    # If use_frequencies is used, this needs to be a dict mapping each frequency to a predict_last_n-value, else an int.
    'predict_last_n': 45,
    # Length of the input sequence
    # If use_frequencies is used, this needs to be a dict mapping each frequency to a seq_length, else an int.
    'seq_length': 365,
    # Length of the forecast sequence 
    'forecast_seq_length': 45,
    # Number of time steps between the forecast issue time and the first forecast
    'forecast_offset': 1,
    # --- Data configurations --------------------------------------------------------------------------
    # which data set to use [camels_us, camels_gb, global, hourly_camels_us, camels_cl, generic]
    'dataset': 'forecast',
    # Path to data set root
    'data_dir': './results/preprocessing/merged',
    # Time series data subdirectory 
    'time_series_data_sub_dir': 'pseudo',
    # Set to True, if train data file should be save to disk. If empty or False, train data is not saved.
    'save_train_data': True,
    # variables to use as time series input (names match the data file column headers)
    'dynamic_inputs': ['rainfall', 'tmean', 'rainfall_forecast', 'tmean_forecast'],
    'forecast_inputs': ['rainfall_forecast', 'tmean_forecast'],
    'hindcast_inputs': ['rainfall', 'tmean'],
    'state_handoff_network': {
        'type': 'fc',
        # define number of neurons per layer in the FC network used as embedding network
        'hiddens': 128,
        # activation function of embedding network
        'activation': 'tanh',
        # dropout applied to embedding network
        'dropout': 0.0
    },
    'target_variables': ['Qd'],
    'clip_targets_to_zero': ['Qd'],
    'static_attributes': [
        'area', 'elev_mean', 'dpsbar', 'sand_perc', 'silt_perc', 'clay_perc',
        'porosity_hypres', 'conductivity_hypres', 'soil_depth_pelletier',
        'frac_snow', 'dwood_perc', 'ewood_perc', 'crop_perc', 'urban_perc',
        'p_mean', 'pet_mean', 'p_seasonality', 'high_prec_freq', 
        'low_prec_freq', 'high_prec_dur', 'low_prec_dur'
    ]
}

# Update configuration with above values
cfg.update_config(cfg_dict)

# Start training 
start_training(cfg=cfg)

# Write `run_dir` to file 
with open(output_file, 'w') as f:
    f.write(str(cfg.run_dir))
