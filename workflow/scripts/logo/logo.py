"""
This script reads in a position weight matrix (PWM), calculates its information content (IC), and generates a logo plot.
"""

import logomaker
from pandas import DataFrame, read_csv
from matplotlib import pyplot as plt

# Snakemake parameters
IP_PWM = snakemake.input[0]  # type: ignore
OP_FIG = snakemake.output[0]  # type: ignore
OP_IC = snakemake.output[1]  # type: ignore
DATASET = snakemake.params[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def read_pwm(filepath: str, names: list[str] = ["A", "C", "G", "T"]) -> DataFrame:
    """Reads input PWM as DataFrame"""
    return read_csv(filepath, header=None, names=names, sep="\t", engine="c")


def calculate_ic(pwm_masked: DataFrame) -> float:
    """Approximates IC from masked PWM"""
    return pwm_masked.max(axis=1).mean()

def calculate_gc(pwm_masked: DataFrame) -> float:
    """Aprroximates GC content from masked PWM"""
    # PWM length
    pwm_length = pwm_masked.shape[0]
    gc_bases = []
    # Loop over PWM
    for row in pwm_masked.itertuples():
        # row info
        base_values = {'A': row.A, 'C': row.C, 'G': row.G, 'T': row.T}
        # Max across row
        row_max = max(base_values.values())
        # Which base is max
        max_base = [base for base, value in base_values.items() if value == row_max][0]
    # If max is g/c
    if max_base == 'G' or max_base == 'C':
        gc_bases.append(1)
    # Return Proportion of GC bases
    return sum(gc_bases) / pwm_length

def main() -> None:
    """
    Main program that reads in a position weight matrix (PWM), calculates its information content (IC),
    generates a logo plot, and saves the IC and the plot.
    """
    # Read pwm
    pwm = read_pwm(IP_PWM)

    # Get index in 1-base
    pwm.index = [i + 1 for i in pwm.index]

    # Mask background odds-ratios
    pwm[pwm < 0] = 0

    # Calculate and store IC approximation
    ic = calculate_ic(pwm)
    
    # Calculate and store GC approximation
    gc = calculate_gc(pwm)

    # Create Logo object
    pwm_logo = logomaker.Logo(pwm, flip_below=False)

    # Style using Logo methods
    pwm_logo.style_spines(visible=False)
    pwm_logo.style_spines(spines=["left", "bottom"], visible=True)
    pwm_logo.style_xticks(rotation=0, fmt="%d", anchor=0)

    # Style using Axes methods
    pwm_logo.ax.set_ylabel("LogOdds", labelpad=-1)
    pwm_logo.ax.xaxis.set_ticks_position("none")
    pwm_logo.ax.xaxis.set_tick_params(pad=-1)
    pwm_logo.ax.set_title(DATASET, color="r")

    # Save IC and GC
    DataFrame({"dataset": [DATASET], "IC": [ic], "GC": [gc]}).to_csv(OP_IC, index=False, sep="\t")

    # Save figure
    plt.savefig(OP_FIG, dpi=100)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
