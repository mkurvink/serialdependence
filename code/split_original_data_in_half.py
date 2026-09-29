from paths import global_path
from pathlib import Path
import pandas as pd

data_n_back_path = global_path / "data_n_back"

def split_by_codenum(X):
    codenum_counts = []

    for path in sorted(Path(data_n_back_path).glob("*.csv")):
        data = pd.read_csv(path, sep=";")

        for (stimulus, codenum), group in data.groupby(["stimulus", "codenum"], sort=False):
            usable_trials = group[f"delta_{X}"].notna().sum()

            codenum_counts.append({
                "source_file": path.stem,
                "stimulus": stimulus,
                "codenum": codenum,
                "usable_trials": usable_trials,
            })

    codenum_counts = pd.DataFrame(codenum_counts)

    assignments = {}

    for stimulus in ["Orientation", "Motion"]:
        stimulus_counts = codenum_counts[
            (codenum_counts["stimulus"] == stimulus)
            & (codenum_counts["usable_trials"] > 0)
        ]

        halfway_point = stimulus_counts["usable_trials"].sum() / 2
        cumulative_trials = 0
        first_half_trials = 0
        second_half_trials = 0

        for _, row in stimulus_counts.iterrows():
            key = (row["source_file"], row["stimulus"], row["codenum"])

            if cumulative_trials < halfway_point:
                assignments[key] = "first_half"
                first_half_trials += row["usable_trials"]
            else:
                assignments[key] = "second_half"
                second_half_trials += row["usable_trials"]

            cumulative_trials += row["usable_trials"]

        print(f"{stimulus}:")
        print(f"  First half:  {first_half_trials} usable trials")
        print(f"  Second half: {second_half_trials} usable trials")
        print(f"  Total:       {first_half_trials + second_half_trials} usable trials")
    
    first_half_path = global_path / "data_n_back_first_half"
    second_half_path = global_path / "data_n_back_second_half"

    first_half_path.mkdir(parents=True, exist_ok=True)
    second_half_path.mkdir(parents=True, exist_ok=True)

    for path in sorted(Path(data_n_back_path).glob("*.csv")):
        data = pd.read_csv(path, sep=";")

        data["split_half"] = [
            assignments.get((path.stem, stimulus, codenum))
            for stimulus, codenum in zip(data["stimulus"], data["codenum"])
        ]

        first_half = data[data["split_half"] == "first_half"].drop(columns=["split_half", "codenum"])
        second_half = data[data["split_half"] == "second_half"].drop(columns=["split_half", "codenum"])

        if len(first_half) > 0:
            first_half.to_csv(first_half_path / path.name, sep=";", index=False)

        if len(second_half) > 0:
            second_half.to_csv(second_half_path / path.name, sep=";", index=False)

    return codenum_counts, assignments