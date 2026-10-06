import subprocess
import sys
from paths import global_path
from n_back_columns import create_n_back_datasets

n_back_path = global_path / "code" / "n_back"

# Step 1: Create the study-specific N-back datasets
create_n_back_datasets()

# Step 2: Create the complete N-back master table in MATLAB
matlab_command = f"addpath('{n_back_path}'); preprocess_create_master_table_n_back"
subprocess.run(["matlab", "-batch", matlab_command], check=True)

# Step 3: Split the resulting master table in Python
subprocess.run(
    [sys.executable, str(n_back_path / "split_master_table.py")],
    check=True,
)

######### Run below only after checking nr. of permutations etc. ########

# # Step 4: Run the MATLAB analysis
# analysis_command = f"addpath('{n_back_path}'); [amplitudes,amp_std_delta,p_values] = estimate_bias_n_back_part_level();"
# subprocess.run(["matlab", "-batch", analysis_command], check=True)