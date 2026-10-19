import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


def plot_timecourses_with_clusters(
    times,
    data,
    clusters=None,
    errors=None,
    labels=None,
    colors=None,
    phases=None,
    events=None,
    title=None,
    xlabel="Time (ms)",
    ylabel="Value",
    figsize=(16, 5.5),
    sigma=0,
    linewidth=2,
    cluster_linewidth=6,
    alpha_error=0.18,
    alpha_phase=0.07,
    cluster_spacing=0.055,
    save_path=None,
    dpi=300,
):
    """
    Plot time-resolved data with significant clusters shown as
    horizontal bars below the time courses.

    Parameters
    ----------
    times : array-like, shape (n_times,)
        Time points, typically in milliseconds.

    data : dict
        Dictionary mapping variable names to 1D time courses.

        Example:
            {
                "Variable A": array_A,
                "Variable B": array_B,
            }

    clusters : dict, optional
        Dictionary mapping variable names to lists of significant
        time intervals (start, end).

        Example:
            {
                "Variable A": [(100, 350), (700, 900)],
                "Variable B": [(400, 600)],
            }

        Only include clusters that have already passed your
        chosen statistical significance threshold.

    errors : dict, optional
        Dictionary mapping variable names to uncertainty arrays,
        such as SEM or confidence interval half-widths.

        Example:
            {
                "Variable A": sem_A,
                "Variable B": sem_B,
            }

    labels : dict, optional
        Display labels for each variable.

    colors : dict, optional
        Colors for each variable.

    phases : list of dict, optional
        Experimental phases to shade.

        Example:
            [
                {
                    "start": 0,
                    "end": 1000,
                    "label": "Phase 1",
                    "color": "tab:blue",
                },
            ]

    events : list of dict, optional
        Vertical event markers.

        Example:
            [
                {
                    "time": 0,
                    "label": "Stimulus onset",
                },
            ]

    sigma : float
        Gaussian smoothing sigma in samples.
        Applied for visualization only.

    save_path : str or Path, optional
        File path for saving the figure.

    Returns
    -------
    fig, ax
        Matplotlib figure and axes.
    """

    times = np.asarray(times, dtype=float)

    if times.ndim != 1:
        raise ValueError("times must be a 1D array.")

    if not data:
        raise ValueError("data must contain at least one time course.")

    clusters = clusters or {}
    errors = errors or {}
    labels = labels or {}
    colors = colors or {}
    phases = phases or []
    events = events or []

    names = list(data.keys())

    # --------------------------------------------------
    # 1. Create figure
    # --------------------------------------------------

    fig, ax = plt.subplots(figsize=figsize)

    # Use default Matplotlib color cycle where needed
    default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    variable_colors = {
        name: colors.get(name, default_colors[i % len(default_colors)])
        for i, name in enumerate(names)
    }

    # --------------------------------------------------
    # 2. Plot experimental phases
    # --------------------------------------------------

    for phase in phases:

        start = phase["start"]
        end = phase["end"]

        phase_color = phase.get("color", "gray")
        phase_label = phase.get("label", "")

        ax.axvspan(
            start,
            end,
            color=phase_color,
            alpha=alpha_phase,
            zorder=0,
        )

        if phase_label:

            midpoint = (start + end) / 2

            ax.text(
                midpoint,
                0.93,
                phase_label,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=10,
                color="0.35",
            )

    # --------------------------------------------------
    # 3. Plot time courses and uncertainty
    # --------------------------------------------------

    all_values = []

    for name in names:

        y = np.asarray(data[name], dtype=float)

        if y.shape != times.shape:
            raise ValueError(
                f"{name}: data shape {y.shape} "
                f"does not match times shape {times.shape}."
            )

        if sigma > 0:
            y_plot = gaussian_filter1d(y, sigma=sigma)
        else:
            y_plot = y

        color = variable_colors[name]

        ax.plot(
            times,
            y_plot,
            label=labels.get(name, name),
            color=color,
            linewidth=linewidth,
            zorder=3,
        )

        all_values.extend(y_plot[np.isfinite(y_plot)])

        if name in errors:

            err = np.asarray(errors[name], dtype=float)

            if err.shape != times.shape:
                raise ValueError(
                    f"{name}: uncertainty shape must match times."
                )

            # Smooth bounds consistently for visualization
            lower = y - err
            upper = y + err

            if sigma > 0:
                lower = gaussian_filter1d(lower, sigma=sigma)
                upper = gaussian_filter1d(upper, sigma=sigma)

            ax.fill_between(
                times,
                lower,
                upper,
                color=color,
                alpha=alpha_error,
                linewidth=0,
                zorder=2,
            )

            all_values.extend(lower[np.isfinite(lower)])
            all_values.extend(upper[np.isfinite(upper)])

    # --------------------------------------------------
    # 4. Determine vertical limits
    # --------------------------------------------------

    if not all_values:
        raise ValueError("No finite values found in data.")

    data_min = min(all_values)
    data_max = max(all_values)

    data_range = data_max - data_min

    if data_range == 0:
        data_range = max(abs(data_max), 1.0)

    # Reserve room for cluster bars below the data
    cluster_names = [
        name for name in names if clusters.get(name)
    ]

    n_cluster_rows = len(cluster_names)

    cluster_top = data_min - 0.10 * data_range

    cluster_positions = {}

    for i, name in enumerate(cluster_names):

        cluster_positions[name] = (
            cluster_top - i * cluster_spacing * data_range
        )

    if cluster_positions:

        bottom = min(cluster_positions.values()) - 0.12 * data_range

    else:

        bottom = data_min - 0.12 * data_range

    top = data_max + 0.15 * data_range

    ax.set_ylim(bottom, top)

    # --------------------------------------------------
    # 5. Draw significant cluster bars
    # --------------------------------------------------

    for name in cluster_names:

        y_position = cluster_positions[name]

        for start, end in clusters[name]:

            ax.hlines(
                y=y_position,
                xmin=start,
                xmax=end,
                color=variable_colors[name],
                linewidth=cluster_linewidth,
                zorder=10,
            )

    # --------------------------------------------------
    # 6. Add event markers
    # --------------------------------------------------

    for event in events:

        event_time = event["time"]
        event_label = event.get("label", "")

        ax.axvline(
            event_time,
            color="0.45",
            linestyle="--",
            linewidth=1,
            alpha=0.75,
            zorder=1,
        )

        if event_label:

            ax.text(
                event_time,
                0.985,
                event_label,
                transform=ax.get_xaxis_transform(),
                rotation=90,
                ha="right",
                va="top",
                fontsize=9,
                color="0.35",
            )

    # --------------------------------------------------
    # 7. Formatting
    # --------------------------------------------------

    ax.set_xlim(times.min(), times.max())

    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)

    if title:
        ax.set_title(title, fontsize=14, pad=12)

    ax.axhline(
        0,
        color="0.5",
        linewidth=0.8,
        alpha=0.5,
        zorder=1,
    )

    ax.tick_params(axis="both", labelsize=11)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ##################### ONLY USE FOR MULTIPLE LINES IN ONE PLOT #############################
    # ax.legend(
    #     loc="upper right",
    #     frameon=False,
    #     fontsize=11,
    # )

    fig.tight_layout()

    # --------------------------------------------------
    # 8. Save
    # --------------------------------------------------

    if save_path is not None:

        fig.savefig(
            save_path,
            dpi=dpi,
            bbox_inches="tight",
        )

    return fig, ax

################################################################################################
################################################################################################
########################  FUNCTION TO USE FOR CLUSTER RECOGNITION  #############################
################################################################################################
################################################################################################

def cluster_recognition(x, clusters):
    p_value_combined = []
    cluster_range_combined = []
    mass_combined = []

    for cluster in clusters:
        for i in range(len(cluster)):
            if (cluster["Lag"][i] == x) & (cluster["P_values"][i] < 0.05):
                p_value = cluster["P_values"][i]
                angles = np.fromstring(cluster["DeltaDegrees"].iloc[i].strip("[]"), sep=",", dtype=int)
                cluster_range = (int(angles[0]), int(angles[-1]))
                mass = cluster["Mass"][i]
                p_value_combined.append(float(p_value))
                cluster_range_combined.append(cluster_range)
                mass_combined.append(float(mass))

    return p_value_combined, cluster_range_combined, mass_combined


################################################################################################
################################################################################################
####################################  EXAMPLE USAGE  ###########################################
################################################################################################
################################################################################################


# # Time points
# times = np.linspace(-500, 2000, 500)

# # Example time courses
# rng = np.random.default_rng(42)

# variable_a = (
#     0.5 * np.exp(-((times - 500) / 300) ** 2)
#     + rng.normal(0, 0.05, len(times))
# )

# variable_b = (
#     -0.4 * np.exp(-((times - 1200) / 250) ** 2)
#     + rng.normal(0, 0.05, len(times))
# )

# # Data to plot
# data = {
#     "Variable A": variable_a,
#     "Variable B": variable_b,
# }

# # Significant clusters (hypothetical)
# clusters = {
#     "Variable A": [
#         (350, 750),
#     ],
#     "Variable B": [
#         (1050, 1400),
#     ],
# }

# # Optional uncertainty
# errors = {
#     "Variable A": np.full(len(times), 0.06),
#     "Variable B": np.full(len(times), 0.06),
# }

# # Experimental phases
# phases = [
#     {
#         "start": 0,
#         "end": 1000,
#         "label": "Phase 1",
#         "color": "tab:blue",
#     },
#     {
#         "start": 1000,
#         "end": 2000,
#         "label": "Phase 2",
#         "color": "tab:orange",
#     },
# ]

# # Event markers
# events = [
#     {
#         "time": 0,
#         "label": "Event 1",
#     },
#     {
#         "time": 1000,
#         "label": "Event 2",
#     },
# ]

# fig, ax = plot_timecourses_with_clusters(
#     times=times,
#     data=data,
#     clusters=clusters,
#     errors=errors,
#     phases=phases,
#     events=events,
#     title="Time-resolved analysis",
#     ylabel="Effect size",
#     sigma=3,                        #### Not sure if this is wanted for the visualization
# )

# plt.show()

# ################################################################################################
# ################################################################################################
# ##############################  IMPLEMENTATION FOR CURRENT DATA  ###############################
# ################################################################################################
# ################################################################################################

from paths import tables_first_half_path, cluster_based_figs_path
import pandas as pd
import numpy as np

# Change in line with data
bin_sizes = [3, 11]
X = 10

for bin_size in bin_sizes:
    print(f"Plotting {bin_size}-degree window figures")
    pos_clusters = pd.read_csv(tables_first_half_path / f"{bin_size}_degree_window" / "positive_clusters.csv", sep = ";")
    np.array(pos_clusters)
    neg_clusters = pd.read_csv(tables_first_half_path / f"{bin_size}_degree_window" / "negative_clusters.csv", sep = ";")
    np.array(neg_clusters)
    save_folder = cluster_based_figs_path / f"{bin_size}_degree_window"
    save_folder.mkdir(parents=True, exist_ok=True)
    for x in range(X):
        p_value_combined, cluster_range_combined, mass_combined = cluster_recognition(x=x+1, clusters = [pos_clusters, neg_clusters])

        group_bias_data = pd.read_csv(tables_first_half_path / f"{bin_size}_degree_window" / f"group_bias_N_{x+1}.csv", sep = ";")
        delta_values = np.array(group_bias_data.loc[:,"DeltaDegrees"])
        group_bias = np.array(group_bias_data.loc[:,"MeanBias"])
        clusters = cluster_range_combined
        save_path = save_folder / f"cluster_based_N_{x+1}"

        events = [{"time": 0}]

        fig, ax = plot_timecourses_with_clusters(
            times = delta_values,
            data={f"N - {x+1}": group_bias},
            clusters={f"N - {x+1}": cluster_range_combined},
            ylabel="Group bias",
            xlabel="Delta values",
            save_path=save_path,
            events=events
        )
        print(f"Finished plotting the figure for N-{x+1}")
