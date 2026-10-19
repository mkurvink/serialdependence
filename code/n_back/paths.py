import pandas as pd
from pathlib import Path

# Paths
global_path = Path("/Users/Kurvi001/Documents/Serial_Dependence/Analysis")
master_path = global_path / "results" / "SD_ma_master_table.csv"
master_n_back_path = global_path / "results" / "SD_ma_master_table_n_back.csv"
data_path = global_path / "data"
results_path = global_path / "results"
code_path = global_path / "code"
figures_first_half_path = results_path / "figures_first_half"
tables_first_half_path = results_path / "tables_first_half"
cluster_based_figs_path = figures_first_half_path / "cluster_based_figs"
