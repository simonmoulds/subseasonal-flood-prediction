# Main entrypoint of the workflow.
# Please follow the best practices:
# https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html,
# in particular regarding the standardized folder structure mentioned there.


rule merge_era5_data:
    """Merge ERA5 data files to make them easier to process"""
    output:
        'results/input/era5_reanalysis_{era5_variable}.nc'
    conda:
        '../envs/conda_environment_base.yml'
    params:
        output_directory='results/input'
    resources:
        time = "00:15:00",
        mem_mb = 8000,
        partition = "short"
    script:
        '../scripts/merge-era5-files.py'


rule extract_era5_field_data:
    """Extract reanalysis climate data for each
    catchment from ERA5 dataset"""
    input:
        'results/input/era5_reanalysis_{era5_variable}.nc'
    output:
        'results/input/ERA5/{station}/{era5_variable}/part-0.parquet'
    params:
        output_directory='results/input/ERA5'
    conda:
        '../envs/conda_environment_r.yml'
    resources:
        time = "02:00:00",
        mem_mb = 8000,
        partition = "short"
    script:
        '../scripts/extract-era5-field-data.R'


rule extract_haduk_field_data:
    """Extract observed climate data for each UK
    catchment from the HadUK dataset"""
    output:
        'results/input/HadUK/{station}/{haduk_variable}/part-0.parquet'
    params:
        output_directory='results/input/HadUK'
    conda:
        '../envs/conda_environment_r.yml'
    resources:
        time = "02:00:00",
        mem_mb = 32000,
        partition = "short"
    script:
        '../scripts/extract-haduk-field-data.R'


rule extract_c3s_field_data:
    """Extract precipitation/temperature for each catchment.

    Only data from model/members specified in the
    configuration file is returned.
    """
    output:
        'results/input/C3S/{station}/{c3s_variable}/part-0.parquet'
    params:
        output_directory='results/input/C3S'
    conda:
        '../envs/conda_environment_r.yml'
    resources:
        time = "01:00:00",
        mem_mb = 32000,
        partition = "short"
    script:
        '../scripts/extract-c3s-field-data.R'










# NOT USED:

# rule extract_hindcast_data:
#     """Compute climate indices for each C3S member"""
#     output:
#         touch('results/done/extract_hindcast_data_{init_month}_{system}.done')
#     conda:
#         'envs/conda_environment_base.yml'
#     resources:
#         time = "12:00:00",
#         mem_mb = 15000,
#         partition = "short"
#     script:
#         'scripts/extract-hindcast-data.py'

# rule extract_eobs_data:
#     """Extract observed climate data for each European
#     catchment from the E-OBS dataset"""
#     input:
#         stations='results/intermediate/station_list_with_grid_coords.parquet',
#     output:
#         'results/intermediate/eobs.parquet'
#     resources:
#         time = "02:00:00",
#         mem_mb = 8000,
#         partition = "short"
#     script:
#         'scripts/extract-eobs-data.R'

# rule extract_reanalysis_data:
#     """Extract reanalysis climate indices from ERA5 dataset"""
#     input:
#         expand('results/intermediate/era5_reanalysis_{clim_variable}.nc', clim_variable=ERA5_VARIABLE)
#     output:
#         'results/intermediate/era5.parquet',
#     conda:
#         'envs/conda_environment_base.yml'
#     resources:
#         time = "02:00:00",
#         mem_mb = 8000,
#         partition = "short"
#     script:
#         'scripts/extract-reanalysis-data.py'
