from paths import results_path, data_path, code_path
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt

def estimate_bias_n_back_part_level():
    tbl = pd.read_csv(results_path / "SD_ma_master_table_n_back_first_half.csv", sep=";")

    # Initializing some variables
    X = 6
    error_variable = "error_ori_deb"
    c = np.sqrt(2) / np.exp(-0.5)
    amp_std_delta = np.full((X, 2), np.nan)
    amp_std_delta_abs = np.full((X, 2), np.nan)
    amplitudes = [None] * X
    fit_p_values = [None] * X
    fdr_p_values = [None] * X
    signed_p_values = np.full(X, np.nan)
    absolute_p_values = np.full(X, np.nan)
    n_permutations = 500    # Number of permutations
    rng = np.random.default_rng(289)  # Random number generator for reproducibility
    # Derivative of Gaussian curve
    def DoG(x, a, peak):
        return x * a * (1/(np.sqrt(2)*peak)) * c * np.exp(-(x/(np.sqrt(2)*peak))**2)
    def fit_dog(x_data,y_data):
        param, covar = scipy.optimize.curve_fit(
            DoG,
            x_data,
            y_data,
            p0=starting_values,
            bounds=(lower_bounds,upper_bounds),
            method="trf",
            max_nfev=2000,
        )

        return param

    # Benjamini-Hochberg procedure for multiple comparisons correction
    def benjamini_hochberg(p_values):
        p_values = np.asarray(p_values,dtype=float)
        adjusted_p_values = np.full(p_values.shape,np.nan)

        valid = ~np.isnan(p_values)
        valid_p_values = p_values[valid]
        m = len(valid_p_values)

        sort_order = np.argsort(valid_p_values)
        sorted_p_values = valid_p_values[sort_order]
        ranks = np.arange(1, m + 1)

        adjusted_sorted_p_values = sorted_p_values * m / ranks
        adjusted_sorted_p_values = np.minimum.accumulate(adjusted_sorted_p_values[::-1])[::-1]
        adjusted_sorted_p_values = np.minimum(adjusted_sorted_p_values,1)

        valid_adjusted_p_values = np.empty(m)
        valid_adjusted_p_values[sort_order] = adjusted_sorted_p_values
        adjusted_p_values[valid] = valid_adjusted_p_values

        return adjusted_p_values

    # Specifications for the curve fitting
    starting_values = [2, 25]  # Initial guesses for the parameters
    lower_bounds = [-20,3]
    upper_bounds = [20,85]

    # Computing the bias for the observed data
    for trial_lat in range(X):
        delta_variable = f"delta_{trial_lat + 1}"
        obs_amplitudes = []
        count_exc = 0
        count_x_data = 0
        studies = np.unique(tbl["studynum"])

        for i in studies:
            tbl_i = tbl.loc[tbl["studynum"] == i]
            tbl_i = tbl_i.loc[(tbl_i[delta_variable] > -90) & (tbl_i[delta_variable] <= 90)]
            experiments = np.unique(tbl_i["expnum"])
            
            for k in experiments:
                tbl_i_k = tbl_i.loc[tbl_i["expnum"] == k]
                cond = np.unique(tbl_i_k["cond"])
                
                for j in cond:
                    tbl_i_k_j = tbl_i_k.loc[tbl_i_k["cond"] == j]
                    obs = np.unique(tbl_i_k_j["obs"])
                    
                    for o in obs:
                        tbl_i_k_j_o = tbl_i_k_j.loc[tbl_i_k_j["obs"] == o]

                        # delta and error values for this specific o
                        signed_delta = tbl_i_k_j_o[delta_variable]
                        signed_error = tbl_i_k_j_o[error_variable]
                        valid = signed_delta.notna() & signed_error.notna()
                        signed_delta = signed_delta.loc[valid].to_numpy()
                        signed_error = signed_error.loc[valid].to_numpy()

                        count_x_data = count_x_data + 1
                        if len(signed_delta) == 0:
                            count_exc = count_exc + 1
                            continue

                        x_data = signed_delta
                        y_data = signed_error

                        if len(np.unique(x_data)) < 5:
                            count_exc = count_exc + 1
                            continue

                        # Fit model
                        param = fit_dog(x_data,y_data)
                        obs_amplitudes.append(param[0])

        print(f"{count_exc} out of {count_x_data} participant-level datasets were excluded for {delta_variable}")

        # Save trial_lat specific amplitudes
        obs_amplitudes = np.asarray(obs_amplitudes)
        amplitudes[trial_lat] = obs_amplitudes
        amp_std_delta[trial_lat,:] = [np.mean(obs_amplitudes),np.std(obs_amplitudes, ddof=1)]
        amp_std_delta_abs[trial_lat,:] = [np.mean(np.abs(obs_amplitudes)), np.std(np.abs(obs_amplitudes), ddof=1)]
        print(f"Mean amplitude for N-{trial_lat + 1}: {amp_std_delta[trial_lat,0]}, SD: {amp_std_delta[trial_lat,1]}")
        print(f"Mean absolute amplitude for N-{trial_lat + 1}: {amp_std_delta_abs[trial_lat,0]}, SD: {amp_std_delta_abs[trial_lat,1]}")

    # # Plotting the mean amplitudes with error bars
    # n_back = np.arange(1, X + 1)

    # fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)

    # axes[0].errorbar(
    #     n_back,
    #     amp_std_delta[:, 0],
    #     yerr=amp_std_delta[:, 1],
    #     marker="o",
    #     capsize=4,
    #     color="black",
    # )
    # axes[0].axhline(0, color="gray", linestyle="--", linewidth=1)
    # axes[0].set_title("Signed amplitudes")
    # axes[0].set_xlabel("Previous trial")
    # axes[0].set_ylabel("Mean amplitude ± SD")

    # axes[1].errorbar(
    #     n_back,
    #     amp_std_delta_abs[:, 0],
    #     yerr=amp_std_delta_abs[:, 1],
    #     marker="o",
    #     capsize=4,
    #     color="black",
    # )
    # axes[1].set_title("Absolute amplitudes")
    # axes[1].set_xlabel("Previous trial")
    # axes[1].set_ylabel("Mean absolute amplitude ± SD")

    # for ax in axes:
    #     ax.set_xticks(n_back)
    #     ax.set_xticklabels([f"N-{i}" for i in n_back])
    #     ax.spines["top"].set_visible(False)
    #     ax.spines["right"].set_visible(False)

    # fig.tight_layout()
    # fig.savefig(results_path / "mean_amplitudes_n_back.png", dpi=300)
    # plt.show()

################################################################################
################## Run only when doing permutations ############################
################################################################################

    # Computing the bias for the observed data
    for trial_lat in range(X):
        delta_variable = f"delta_{trial_lat + 1}"
        permuted_obs_amplitudes = []
        observed_amplitudes = np.asarray(amplitudes[trial_lat])
        n_observed_fits = len(observed_amplitudes)

        fit_perm_higher_counts = np.zeros(n_observed_fits, dtype=int)
        permuted_signed_means = np.full(n_permutations, np.nan)
        permuted_absolute_means = np.full(n_permutations, np.nan)

        for permutation in range(n_permutations):
            obs_amplitudes = []
            studies = np.unique(tbl["studynum"])

            for i in studies:
                tbl_i = tbl.loc[tbl["studynum"] == i]
                tbl_i = tbl_i.loc[(tbl_i[delta_variable] > -90) & (tbl_i[delta_variable] <= 90)]
                experiments = np.unique(tbl_i["expnum"])
                
                for k in experiments:
                    tbl_i_k = tbl_i.loc[tbl_i["expnum"] == k]
                    cond = np.unique(tbl_i_k["cond"])
                    
                    for j in cond:
                        tbl_i_k_j = tbl_i_k.loc[tbl_i_k["cond"] == j]
                        obs = np.unique(tbl_i_k_j["obs"])
                        
                        for o in obs:
                            tbl_i_k_j_o = tbl_i_k_j.loc[tbl_i_k_j["obs"] == o]

                            # delta and error values for this specific o
                            signed_delta = tbl_i_k_j_o[delta_variable]
                            signed_error = tbl_i_k_j_o[error_variable]
                            valid = signed_delta.notna() & signed_error.notna()
                            signed_delta = signed_delta.loc[valid].to_numpy()
                            signed_error = signed_error.loc[valid].to_numpy()

                            if len(signed_delta) == 0:
                                continue

                            if len(np.unique(x_data)) < 5:
                                continue

                            shuffled_delta = rng.permutation(signed_delta)

                            x_data = shuffled_delta
                            y_data = signed_error

                            # Fit model
                            param = fit_dog(x_data,y_data)
                            obs_amplitudes.append(param[0])

            permuted_obs_amplitudes = np.asarray(obs_amplitudes)
            if len(permuted_obs_amplitudes) != n_observed_fits:
                raise RuntimeError(f"Number of fits differs during permutation {permutation + 1} "f"for N-{trial_lat + 1}")

            fit_perm_higher_counts += np.abs(permuted_obs_amplitudes) >= np.abs(observed_amplitudes)   # Adds a 1 if the permuted amplitude is greater than or equal to the observed amplitude, resulting in f.e. [0, 534, 82, 25, 10, 5245,...]
            permuted_signed_means[permutation] = np.mean(permuted_obs_amplitudes)
            permuted_absolute_means[permutation] = np.mean(np.abs(permuted_obs_amplitudes))

            if (permutation + 1) % (n_permutations / 10) == 0:
                print(f"N-{trial_lat + 1}: {permutation + 1}/{n_permutations} permutations completed")

        # Calculate p-value per trial_lat 
        fit_p_values[trial_lat] = (1 + fit_perm_higher_counts) / (n_permutations + 1)
        fdr_p_values[trial_lat] = benjamini_hochberg(fit_p_values[trial_lat])
        significant_fdr = fdr_p_values[trial_lat] <= 0.05

        positive = observed_amplitudes > 0
        negative = observed_amplitudes < 0
        positive_count = np.sum(positive)
        negative_count = np.sum(negative)
        significant_positive_count = np.sum(significant_fdr & positive)
        significant_negative_count = np.sum(significant_fdr & negative)

        print(f"N-{trial_lat + 1} fit-level results:\n")
        print(f"Positive amplitudes: {positive_count}/{n_observed_fits}\n")
        print(f"Negative amplitudes: {negative_count}/{n_observed_fits}\n")
        print(f"Significant positive amplitudes: {significant_positive_count}/{n_observed_fits}\n")
        print(f"Significant negative amplitudes: {significant_negative_count}/{n_observed_fits}\n")

        observed_signed_mean = np.mean(observed_amplitudes)
        observed_absolute_mean = np.mean(np.abs(observed_amplitudes))
        signed_p_values[trial_lat] = (1 + np.sum(np.abs(permuted_signed_means) >= np.abs(observed_signed_mean))) / (n_permutations + 1)
        absolute_p_values[trial_lat] = (1 + np.sum(permuted_absolute_means >= observed_absolute_mean)) / (n_permutations + 1)

        print(f"Signed mean amplitude p-value for N-{trial_lat + 1}: {signed_p_values[trial_lat]}")
        print(f"Absolute mean amplitude p-value for N-{trial_lat + 1}: {absolute_p_values[trial_lat]}")

    fit_p_table = pd.DataFrame({f"delta_{lag+1}": pd.Series(fit_p_values[lag]) for lag in range(X)})
    fit_p_table.to_csv(results_path / "fit_p_values_n_back_first_half.csv", sep = ";", index = False)

    return (amplitudes, amp_std_delta, fit_p_values, signed_p_values, absolute_p_values)

if __name__ == "__main__":
    estimate_bias_n_back_part_level()


