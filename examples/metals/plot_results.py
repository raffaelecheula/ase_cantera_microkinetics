# -------------------------------------------------------------------------------------
# IMPORTS
# -------------------------------------------------------------------------------------

import yaml
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------------------------

def main():

    # Parameters.
    surface_name = "Rh100"
    results_yaml = f"results_{surface_name}.yaml"
    
    # Read input yaml files.
    with open(results_yaml, "r") as fileobj:
        results_all = yaml.safe_load(fileobj)
    
    # Plot log of TOF.
    log_TOF = np.log(results_all["TOF"])
    mean_log_TOF = np.mean(log_TOF)
    std_log_TOF = np.std(log_TOF)
    print(f"log(TOF [1/s]) = {mean_log_TOF:.2f} ± {std_log_TOF:.2f}")
    fig, ax = plt.subplots()
    violins = ax.violinplot(log_TOF, showmeans=False, showextrema=False)
    for body in violins["bodies"]:
        body.set_facecolor("crimson")
        body.set_edgecolor("black")
        body.set_alpha(0.7)
    ax.errorbar(1, mean_log_TOF, yerr=std_log_TOF, fmt="o", capsize=5, color="black")
    ax.set_ylabel("log(TOF [1/s])")
    ax.set_xticks([])
    ax.set_ylim([-30, +10])
    ax.set_xlim([0.5, 1.5])
    ax.set_title(f"log(TOF [1/s]) = {mean_log_TOF:.2f} ± {std_log_TOF:.2f}")
    # Save plot.
    plt.tight_layout()
    plt.savefig(f"plot_log_TOF_{surface_name}.png", dpi=300)

    # Plot coverages.
    fig, ax = plt.subplots()
    keys_theta = [key for key in results_all.keys() if key.startswith("theta")]
    results_theta = [results_all[key] for key in keys_theta]
    violins = ax.violinplot(results_theta, showmeans=False, showextrema=False)
    for body in violins["bodies"]:
        body.set_facecolor("crimson")
        body.set_edgecolor("black")
        body.set_alpha(0.7)
    ax.set_ylabel("Coverages [-]")
    ax.set_xticks(range(1, len(keys_theta) + 1))
    names_theta = [cantera_to_labels(key.replace("theta", "")) for key in keys_theta]
    ax.set_xticklabels(names_theta, rotation=90)
    ax.set_ylim([0, 1])
    # Save plot.
    plt.tight_layout()
    plt.savefig(f"plot_coverages_{surface_name}.png", dpi=300)

    # Plot DRCs.
    fig, ax = plt.subplots()
    keys_DRC = [key for key in results_all.keys() if key.startswith("DRC")]
    results_DRC = [results_all[key] for key in keys_DRC]
    results_DRC = [
        [value for value in result if 0 < value < 1] for result in results_DRC
    ]
    violins = ax.violinplot(results_DRC, showmeans=False, showextrema=False)
    for body in violins["bodies"]:
        body.set_facecolor("crimson")
        body.set_edgecolor("black")
        body.set_alpha(0.7)
    ax.set_ylabel("Degree of Rate Control [-]")
    ax.set_xticks(range(1, len(keys_DRC) + 1))
    names_DRC = [cantera_to_labels(key.replace("DRC", "")) for key in keys_DRC]
    ax.set_xticklabels(names_DRC, rotation=90)
    ax.set_ylim([0, 1])
    # Save plot.
    plt.tight_layout()
    plt.savefig(f"plot_DRCs_{surface_name}.png", dpi=300)

# -------------------------------------------------------------------------------------
# CANTERA TO LABELS
# -------------------------------------------------------------------------------------

def cantera_to_labels(
    label: str,
) -> str:
    """
    Convert Cantera key to label for plotting.
    """
    replacements = {
        "O2": "O$_2$",
        "H2": "H$_2$",
        "(Rh)": "*",
        "(Rh,Rh)": "**",
        "<=>": "→",
        " ": "",
    }
    for key, value in replacements.items():
        label = label.replace(key, value)
    return label

# -------------------------------------------------------------------------------------
# IF NAME MAIN
# -------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()

# -------------------------------------------------------------------------------------
# END
# -------------------------------------------------------------------------------------
