import pandas as pd
import numpy as np
from pathlib import Path
from paths import global_path, master_path, data_path

# Paths
df = pd.read_csv(master_path, sep=",")  # The comma/ semicolon separator may sometimes change somehow; so for any random error change the sep argument
output_path = global_path / "data_n_back"
output_path.mkdir(parents=True, exist_ok=True)

def create_n_back_datasets(output_path = output_path, data_path = data_path, df = df):
    source_dfs = {}
    source_master_dfs = {}
    codenum_summaries = []

    for path in sorted(Path(data_path).glob("*.csv")):
        name = path.stem
        source_dfs[name] = pd.read_csv(path, sep=";")
        original_columns = source_dfs[name].columns.tolist()
        study_name = source_dfs[name]["study"].iloc[0]
        print(name)

        # Reproduce the master-table preprocessing for missing conditions or blocks
        if "cond" not in source_dfs[name].columns:
            source_dfs[name]["cond"] = 1
        if "block" not in source_dfs[name].columns:
            source_dfs[name]["block"] = 1

        source_master_dfs[name] = (
            df[df["study"] == study_name]
            .copy()
            .reset_index(drop=True)
        )

        # Preserve the row order
        source_dfs[name]["_row"] = np.arange(len(source_dfs[name]))
        source_dfs[name] = (source_dfs[name].sort_values(["expnum", "cond", "obs", "_row"], kind="stable").reset_index(drop=True))

        # Check that the source and master subsets contain the same total number of trials.
        assert len(source_dfs[name]) == len(source_master_dfs[name])

        # Condition, block, and stimulus should agree row by row.
        assert np.array_equal(source_dfs[name]["cond"].to_numpy(), source_master_dfs[name]["cond"].to_numpy())
        assert np.array_equal(source_dfs[name]["block"].to_numpy(), source_master_dfs[name]["block"].to_numpy())
        assert np.array_equal(source_dfs[name]["stimulus"].to_numpy(), source_master_dfs[name]["stimulus"].to_numpy())

        # Find the previous retained trial
        # Variables that define the participant-level experimental context
        context_cols = ["study", "obs", "cond", "expnum"]

        # Within each context, determine whether the block number has changed
        block_changed = (source_dfs[name].groupby(context_cols, sort=False)["block"].transform(lambda block: block.ne(block.shift())))

        # Give every continuous occurrence of a block a new number
        source_dfs[name]["block_unique"] = (
            block_changed.groupby(
                [source_dfs[name][column] for column in context_cols],
                sort=False,
            ).cumsum())

        # Use study, original participant number, block, and condition as grouping variables.
        group_cols = ["study", "obs", "cond", "block_unique", "expnum"]
        g = source_dfs[name].groupby(group_cols, sort=False)

        # Relevant for applying orientation wrapping rule
        orientation = source_dfs[name]["stimulus"] == "Orientation"
        motion = source_dfs[name]["stimulus"] == "Motion"

        # Def for rounding delta's theta's, ensuring 0.5 is rounded up, and -0.5 is rounded down, as in MATLAB's round function
        def matlab_round(values):
            # Remove tiny floating-point errors from subtraction
            cleaned_values = values.round(10)

            return np.sign(cleaned_values) * np.floor(np.abs(cleaned_values) + 0.5)

        # Choose nr. of N-back trials
        X = 6  # Change this value to select the desired N-back trial count      #####################################                                       ##############################################################
        for i in range(1, X + 1):
            source_dfs[name][f"previous_theta_{i}"] = g["theta"].shift(i)

            # Calculate delta before rounding
            source_dfs[name][f"delta_{i}"] = (source_dfs[name][f"previous_theta_{i}"] - source_dfs[name]["theta"])

            # Apply the orientation wrapping rule
            source_dfs[name].loc[orientation & (source_dfs[name][f"delta_{i}"] > 90), f"delta_{i}"] -= 180
            source_dfs[name].loc[orientation & (source_dfs[name][f"delta_{i}"] < -90), f"delta_{i}"] += 180

            # Apply the equivalent motion wrapping rule
            source_dfs[name].loc[motion & (source_dfs[name][f"delta_{i}"] > 180), f"delta_{i}"] -= 360
            source_dfs[name].loc[motion & (source_dfs[name][f"delta_{i}"] < -180), f"delta_{i}"] += 360

            # Rounding delta's and theta's
            source_dfs[name][f"rounded_delta_{i}"] = matlab_round(source_dfs[name][f"delta_{i}"])
            source_dfs[name][f"rounded_previous_theta_{i}"] = matlab_round(source_dfs[name][f"previous_theta_{i}"])
            if "rounded_theta" not in source_dfs[name]:
                source_dfs[name]["rounded_theta"] = matlab_round(source_dfs[name]["theta"])

            # Match the master table's -90/+90 convention
            source_dfs[name].loc[(source_dfs[name][f"rounded_delta_{i}"] == -90), f"rounded_delta_{i}"] = 90

        # Add the corresponding master-table information
        source_dfs[name]["codenum"] = source_master_dfs[name]["codenum"]
        source_dfs[name]["obsid"] = source_master_dfs[name]["obsid"]
        source_dfs[name]["master_delta"] = source_master_dfs[name]["delta"]
        source_dfs[name]["master_theta"] = source_master_dfs[name]["theta"]

        # Compare the newly calculated delta with the master delta
        source_dfs[name]["testable"] = (source_dfs[name]["previous_theta_1"].notna() & source_dfs[name]["theta"].notna() & source_dfs[name]["master_delta"].notna())
        exact_match = (source_dfs[name]["rounded_delta_1"] == source_dfs[name]["master_delta"])
        equivalent_180 = (source_dfs[name]["rounded_delta_1"].abs().eq(180) & source_dfs[name]["master_delta"].abs().eq(180))
        source_dfs[name]["match"] = (source_dfs[name]["testable"] & (exact_match | equivalent_180))

        # To add new trial count per unique block, to use in the future
        source_dfs[name]["trial_in_block"] = (g.cumcount() + 1)

        codenum_summary = (
            source_dfs[name]
            .groupby(["study", "codenum"], sort=False)
            .agg(
                total_rows=("master_delta", "size"),
                comparisons=("testable", "sum"),
                matches=("match", "sum"),
            )
            .reset_index()
        )
        codenum_summary["mismatches"] = (codenum_summary["comparisons"] - codenum_summary["matches"])
        codenum_summary["passes_delta_check"] = (codenum_summary["mismatches"] == 0)

        codenum_summary["source_file"] = name
        codenum_summaries.append(codenum_summary)

        # Renaming columns to my liking
        pass_by_codenum = codenum_summary.set_index("codenum")["passes_delta_check"]
        source_dfs[name]["passes_delta_check"] = source_dfs[name]["codenum"].map(pass_by_codenum)
        assert source_dfs[name]["passes_delta_check"].notna().all()

        delta_cols = []
        for i in range(1, X + 1):
            delta_col = f"delta_{i}"
            source_dfs[name][delta_col] = source_dfs[name][f"rounded_delta_{i}"]
            delta_cols.append(delta_col)

        failing_codenum = ~source_dfs[name]["passes_delta_check"]
        too_early_in_block = source_dfs[name]["trial_in_block"] <= X
        missing_any_delta = source_dfs[name][delta_cols].isna().any(axis=1)

        set_to_nan = failing_codenum | too_early_in_block | missing_any_delta
        source_dfs[name].loc[set_to_nan, delta_cols] = np.nan

        columns_to_save = original_columns + ["codenum"] + delta_cols + ["trial_in_block", "block_unique"]

        output_df = source_dfs[name].sort_values("_row").reset_index(drop=True)
        output_df = output_df[columns_to_save].copy()
        output_df.to_csv(output_path / path.name, sep=";", index=False)

    all_codenum_summary = pd.concat(codenum_summaries, ignore_index=True)
    summary_path = global_path / "results" / "delta_mismatch_check_by_codenum.csv"
    all_codenum_summary.to_csv(summary_path, index=False)
    print(f"Saved codenum summary to {summary_path}")

    return all_codenum_summary, source_dfs

