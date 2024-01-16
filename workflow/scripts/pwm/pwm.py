#!/usr/bin/python

from pathlib import Path
from Bio import motifs
from pandas import read_csv
from numpy import int64


# Snakemake parameters
IP_MATRIX = snakemake.input[0]  # type: ignore
OP_MATRIX = snakemake.output[0]  # type: ignore

###
# Functions
###


def input_format(filepath: str) -> str:
    """Return format of input motif matrix"""
    return Path(filepath).suffix.split(".")[-1]


def pwm_to_IntLogOdds(filepath: str, outpath: str) -> None:
    """Converts PWM to integer log-odds format"""
    # Read PWM, values are assumed to be log-odds
    pwm = read_csv(filepath, header=None, sep="\t")
    # Round to integer, transpose for PWMScan format
    pwm = pwm.round().astype(int64).T
    # Write out
    pwm.to_csv(outpath, header=None, sep="\t", index=False)


def jaspar_to_IntLogOdds(filepath: str, outpath: str):
    """Converts JASPAR to integer log-odds format"""
    # Read PFM in jaspar format
    motif = motifs.read(filepath, "jaspar")

    # Calculate pseudocounts
    motif.pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)  # type: ignore

    # Caclulate PWM
    pwm = list(map(list, zip(*[motif.pssm[nt] for nt in "ACGT"])))  # type: ignore

    # Write out PWM
    with open(outpath, mode="w", encoding="utf-8") as op:
        for line in pwm:
            line = "\t".join(
                ["{:7d}".format(round(j * 100)) for j in line]
            )  # type: ignore
            op.write(f"{line}\n")


def main() -> None:
    """"""
    # Determine input format
    ip_format = input_format(IP_MATRIX)

    # Run program
    if ip_format == "pwm":
        pwm_to_IntLogOdds(IP_MATRIX, OP_MATRIX)
    if ip_format == "jaspar":
        jaspar_to_IntLogOdds(IP_MATRIX, OP_MATRIX)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
