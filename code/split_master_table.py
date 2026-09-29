import pandas as pd
from paths import results_path

df = pd.read_csv(results_path / "SD_ma_master_table_n_back.csv", sep=";")
usable_codenums = (df.groupby("codenum")["delta_4"].apply(lambda values: values.notna().any()))
usable_codenums = usable_codenums[usable_codenums].index
df = df[df["codenum"].isin(usable_codenums)].copy()

first_half_parts = []
second_half_parts = []

for stimulus in df["stimulus"].unique():
    stimulus_data = df[df["stimulus"] == stimulus]

    halfway_point = len(stimulus_data) / 2
    cumulative_trials = 0

    for obsid in stimulus_data["obsid"].unique():
        participant_data = stimulus_data[stimulus_data["obsid"] == obsid]

        if cumulative_trials < halfway_point:
            first_half_parts.append(participant_data)
            cumulative_trials += len(participant_data)
        else:
            second_half_parts.append(participant_data)

    first_trials = sum(len(part) for part in first_half_parts if part["stimulus"].iloc[0] == stimulus)
    second_trials = sum(len(part) for part in second_half_parts if part["stimulus"].iloc[0] == stimulus)

    print(f"{stimulus}:")
    print(f"  First half:  {first_trials} trials")
    print(f"  Second half: {second_trials} trials")

first_half = pd.concat(first_half_parts, ignore_index=True)
second_half = pd.concat(second_half_parts, ignore_index=True)

first_half_path = results_path / "SD_ma_master_table_n_back_first_half.csv"
second_half_path = results_path / "SD_ma_master_table_n_back_second_half.csv"

first_half.to_csv(first_half_path, sep=";", index=False)
second_half.to_csv(second_half_path, sep=";", index=False)

print(f"\nSaved first half: {first_half_path}")
print(f"Saved second half: {second_half_path}")
print(f"First-half rows: {len(first_half)}")
print(f"Second-half rows: {len(second_half)}")