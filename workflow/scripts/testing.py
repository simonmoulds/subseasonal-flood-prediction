import pickle

from pathlib import Path

from neuralhydrology.evaluation import get_tester
from neuralhydrology.utils.config import Config

# Snakemake stuff
input_filename = Path(snakemake.input[0])
output_filename = Path(snakemake.output[0])
member = str(snakemake.wildcards['member'])

# # TESTING
# input_filename = Path('results/training/latest_run_dir.txt')
# output_filename = Path('results/testing/test.p')
# member = '1'

with open(input_filename, 'r') as f:
    run_dir = Path(f.read().strip())

# https://github.com/neuralhydrology/neuralhydrology/discussions/107#discussioncomment-6198244
def start_prediction(cfg: Config, run_dir: Path, time_series_sub_dir: str, output_filename: Path, epoch: int = None):
    """Start prediction for a list of basins over a defined date range using a trained network. 
    
    Useful for generation predictions at ungaged locations where no observations are available.  Does not assume observations are available so no metrics are saved

    Parameters
    ----------
    cfg : Config
        The run configuration, read from the run directory.
    run_dir : Path
        Path to the run directory.
    time_series_sub_dir : str 
        Subdirectory containing forcing directory
    output_filename : Path
        Path to output file for prediction results
    epoch : int, optional
        Define a specific epoch to evaluate. By default, the weights of the last epoch are used.
    """

    cfg.update_config({'time_series_data_sub_dir': time_series_sub_dir})
    tester = get_tester(cfg=cfg, run_dir=run_dir, period = 'test', init_model=True)
    pred_results = tester.evaluate(epoch=epoch, save_results=False)

    # file_name = output_dir / time_series_sub_dir / f"{Path(cfg.test_basin_file).stem}_prediction_results.p"
    # file_name.parent.mkdir(parents=True, exist_ok=True)
    with output_filename.open("wb") as fp:
        pickle.dump(pred_results, fp)

cfg = Config(run_dir / "config.yml")
start_prediction(cfg=cfg, run_dir=run_dir, time_series_sub_dir=member, output_filename=output_filename, epoch=30)

# # TODO this should be the next stage of the workflow:
# # TESTING
# # This shows how to extract data for one station / member combination, retrieve corresponding EFAS data, and join together
# import pandas as pd
# import pickle 
# with open('results/testing/prediction/prediction_results_4.p', 'rb') as f:
#     data = pickle.load(f)

# ids = list(data.keys())
# id = ids[0] # TESTING ONLY
# xr = data[id]['1D']['xr']
# df = xr.to_dataframe()
# df = df.dropna()
# df = df.reset_index()
# last_time = df['date']
# time_step = df['time_step']
# time = [last_time[i] + pd.Timedelta(days=int(time_step[i] + 1)) for i in range(len(last_time))]
# max_time_step = time_step.min()
# init_time = [last_time[i] + pd.Timedelta(days=int(max_time_step)) for i in range(len(last_time))]
# lead_time = [(time[i] - init_time[i]).days for i in range(len(time))]
# df['member'] = member
# df['time'] = time 
# df['init_time'] = init_time 
# df['lead_time'] = lead_time 
# df = df[['member', 'init_time', 'lead_time', 'time', 'Qd_obs', 'Qd_sim']]

# # Also compare with EFAS 
# efas_df = pd.read_parquet(Path('results/preprocessing/EFAS') / f'{id}_bc.parquet')
# efas_df['init_time'] = pd.to_datetime(efas_df['init_time'])
# efas_df['time'] = pd.to_datetime(efas_df['time'])
# efas_df = efas_df[['member', 'init_time', 'lead_time', 'time', 'value_bc']]
# efas_df = efas_df.rename({'value_bc': 'Qd_efas'}, axis=1)

# # Join simulated data with EFAS data
# newdf = df.merge(efas_df, how='left', on=['member', 'init_time', 'lead_time', 'time'])

# NOT USED:
# Vector = list[str]
# def start_ensemble_prediction(cfg: Config, run_dir: Path, time_series_sub_dirs: Vector, output_dir: Path, epoch: int = None): 
#     pred_results_dict = {}
#     for sub_dir in time_series_sub_dirs:
#         cfg.update_config({'time_series_data_sub_dir': sub_dir})
#         tester = get_tester(cfg=cfg, run_dir=run_dir, period = 'test', init_model=True)
#         pred_results = tester.evaluate(epoch=epoch, save_results=False)
#         pred_results_dict[sub_dir] = pred_results
#     file_name = output_dir / f"{Path(cfg.test_basin_file).stem}_prediction_results.p"
#     file_name.parent.mkdir(parents=True, exist_ok=True)
#     with file_name.open("wb") as fp:
#         pickle.dump(pred_results_dict, fp)
# members = [str(i) for i in range(0, 24)]
# start_ensemble_prediction(cfg, run_dir=run_dir, time_series_sub_dirs=members, output_dir=output_dir, epoch=30)
# # Test data: 
# fname = output_dir / "basin_list_prediction_results.p"
# with open(fname, 'rb') as file: 
#     x = pickle.load(file)
# station = '10002'
# # Plotting
# plt.figure(figsize=(10, 6))
# for key in x:
#     key = '10'
#     xr = x[key][station]['1D']['xr']
#     xr = xr.isel(date=0)
#     df = xr.to_dataframe()
#     dates = df['date'].values 
#     dates = [dates[i] + pd.Timedelta(days=int(df.index[i] + 1)) for i in range(df.shape[0])]
#     df['time'] = dates 
#     plt.plot(df['time'], df['Qd_sim'], label=None, color='grey')
# plt.plot(df['time'], df['Qd_obs'], label='Obs', color='b')
# y = pd.read_parquet("data/EFAS/20081001.parquet")
# y = y[y['id'] == station]
# y = y[y['member'] == '10']
# y = y.iloc[0:45]
# plt.plot(y['time'], y['value'], label=None, color='g')
# plt.xlabel('Time')
# plt.ylabel('Values')
# plt.title('Time Series Plot of Variable 1 and Variable 2')
# plt.legend()
# plt.grid(True)
# # plt.show()
# plt.savefig('test.png')

