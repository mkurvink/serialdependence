from paths import global_path, results_path
import pandas as pd
import numpy as np
from scipy.stats import pearsonr

amplitude_path = results_path / "amplitudes_n_back_first_half.csv"
df = pd.read_csv(amplitude_path, sep=";", skiprows=1, header=None)
row_length = len(df)
column_length = len(df.columns)
df_binary = df.copy()
positive_N_1 = pd.DataFrame(columns=df.columns)

for i in range(row_length):
    for j in range(column_length):
        if df.loc[i,j] > 0:
            df_binary.loc[i,j] = 1
        else:
            df_binary.loc[i,j] = 0

count_n12 = count_n13 = count_n14 = count_n23 = count_n24 = count_n34 = 0
for i in range(row_length):
    if (df_binary.loc[i,0] == 1.0) & (df_binary.loc[i,1] == 1.0):
        count_n12 += 1
    if (df_binary.loc[i,0] == 1.0) & (df_binary.loc[i,2] == 1.0):
        count_n13 += 1
    if (df_binary.loc[i,0] == 1.0) & (df_binary.loc[i,3] == 1.0):
        count_n14 += 1
    if (df_binary.loc[i,1] == 1.0) & (df_binary.loc[i,2] == 1.0):
        count_n23 += 1
    if (df_binary.loc[i,1] == 1.0) & (df_binary.loc[i,3] == 1.0):
        count_n24 += 1
    if (df_binary.loc[i,2] == 1.0) & (df_binary.loc[i,3] == 1.0):
        count_n34 += 1

print(f"Positive overlap between N - 1 and N - 2: {count_n12}")
print(f"Positive overlap between N - 1 and N - 3: {count_n13}")
print(f"Positive overlap between N - 1 and N - 4: {count_n14}")
print(f"Positive overlap between N - 2 and N - 3: {count_n23}")
print(f"Positive overlap between N - 2 and N - 4: {count_n24}")
print(f"Positive overlap between N - 3 and N - 4: {count_n34}\n")

neg_count_n12 = neg_count_n13 = neg_count_n14 = neg_count_n23 = neg_count_n24 = neg_count_n34 = 0
for i in range(row_length):
    if (df_binary.loc[i,0] == 0.0) & (df_binary.loc[i,1] == 0.0):
        neg_count_n12 += 1
    if (df_binary.loc[i,0] == 0.0) & (df_binary.loc[i,2] == 0.0):
        neg_count_n13 += 1
    if (df_binary.loc[i,0] == 0.0) & (df_binary.loc[i,3] == 0.0):
        neg_count_n14 += 1
    if (df_binary.loc[i,1] == 0.0) & (df_binary.loc[i,2] == 0.0):
        neg_count_n23 += 1
    if (df_binary.loc[i,1] == 0.0) & (df_binary.loc[i,3] == 0.0):
        neg_count_n24 += 1
    if (df_binary.loc[i,2] == 0.0) & (df_binary.loc[i,3] == 0.0):
        neg_count_n34 += 1

print(f"Negative overlap between N - 1 and N - 2: {neg_count_n12}")
print(f"Negative overlap between N - 1 and N - 3: {neg_count_n13}")
print(f"Negative overlap between N - 1 and N - 4: {neg_count_n14}")
print(f"Negative overlap between N - 2 and N - 3: {neg_count_n23}")
print(f"Negative overlap between N - 2 and N - 4: {neg_count_n24}")
print(f"Negative overlap between N - 3 and N - 4: {neg_count_n34}\n")

count_n123 = count_n1234 = count_n234 = 0
for i in range(row_length):
    if (df_binary.loc[i,0] == 1.0) & (df_binary.loc[i,1] == 1.0) & (df_binary.loc[i,2] == 1.0):
        count_n123 += 1
    if (df_binary.loc[i,0] == 1.0) & (df_binary.loc[i,1] == 1.0) & (df_binary.loc[i,2] == 1.0) & (df_binary.loc[i,3] == 1.0):
        count_n1234 += 1
    if (df_binary.loc[i,1] == 1.0) & (df_binary.loc[i,2] == 1.0) & (df_binary.loc[i,3] == 1.0):
        count_n234 += 1

print(f"Positive overlap between N - 1, N - 2, and N - 3: {count_n123}")
print(f"Positive overlap between N - 1, N - 2, N - 3, and N - 4: {count_n1234}")
print(f"Positive overlap between N - 2, N - 3, and N - 4: {count_n234}\n")

neg_count_n123 = neg_count_n1234 = neg_count_n234 = 0
for i in range(row_length):
    if (df_binary.loc[i,0] == 0.0) & (df_binary.loc[i,1] == 0.0) & (df_binary.loc[i,2] == 0.0):
        neg_count_n123 += 1
    if (df_binary.loc[i,0] == 0.0) & (df_binary.loc[i,1] == 0.0) & (df_binary.loc[i,2] == 0.0) & (df_binary.loc[i,3] == 0.0):
        neg_count_n1234 += 1
    if (df_binary.loc[i,1] == 0.0) & (df_binary.loc[i,2] == 0.0) & (df_binary.loc[i,3] == 0.0):
        neg_count_n234 += 1

print(f"Negative overlap between N - 1, N - 2, and N - 3: {neg_count_n123}")
print(f"Negative overlap between N - 1, N - 2, N - 3, and N - 4: {neg_count_n1234}")
print(f"Negative overlap between N - 2, N - 3, and N - 4: {neg_count_n234}\n")


print(f"Correlation between N - 1 and N - 2: {np.corrcoef(df_binary.iloc[:, 0], df_binary.iloc[:, 1])[0, 1]}")
print(f"Correlation between N - 1 and N - 3: {np.corrcoef(df_binary.iloc[:, 0], df_binary.iloc[:, 2])[0, 1]}")
print(f"Correlation between N - 1 and N - 4: {np.corrcoef(df_binary.iloc[:, 0], df_binary.iloc[:, 3])[0, 1]}")
print(f"Correlation between N - 2 and N - 3: {np.corrcoef(df_binary.iloc[:, 1], df_binary.iloc[:, 2])[0, 1]}")
print(f"Correlation between N - 2 and N - 4: {np.corrcoef(df_binary.iloc[:, 1], df_binary.iloc[:, 3])[0, 1]}")
print(f"Correlation between N - 3 and N - 4: {np.corrcoef(df_binary.iloc[:, 2], df_binary.iloc[:, 3])[0, 1]}\n")

r_n12, p_n12 = pearsonr(df_binary.iloc[:, 0], df_binary.iloc[:, 1])
r_n13, p_n13 = pearsonr(df_binary.iloc[:, 0], df_binary.iloc[:, 2])
r_n14, p_n14 = pearsonr(df_binary.iloc[:, 0], df_binary.iloc[:, 3])
r_n23, p_n23 = pearsonr(df_binary.iloc[:, 1], df_binary.iloc[:, 2])
r_n24, p_n24 = pearsonr(df_binary.iloc[:, 1], df_binary.iloc[:, 3])
r_n34, p_n34 = pearsonr(df_binary.iloc[:, 2], df_binary.iloc[:, 3])

print(f"Correlation between N-1 and N-2: r = {r_n12:.3f}, p = {p_n12:.4f}")
print(f"Correlation between N-1 and N-3: r = {r_n13:.3f}, p = {p_n13:.4f}")
print(f"Correlation between N-1 and N-4: r = {r_n14:.3f}, p = {p_n14:.4f}")
print(f"Correlation between N-2 and N-3: r = {r_n23:.3f}, p = {p_n23:.4f}")
print(f"Correlation between N-2 and N-4: r = {r_n24:.3f}, p = {p_n24:.4f}")
print(f"Correlation between N-3 and N-4: r = {r_n34:.3f}, p = {p_n34:.4f}")