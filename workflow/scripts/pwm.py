#!/usr/bin/python

from Bio import motifs

# Snakemake parameters
INPUT_PFM = snakemake.input[0]  # type: ignore
OUTPUT_PWM = snakemake.output[0] # type: ignore

###
# Functions
###

def main():
    """D"""
    print("BARKK")
    with open(INPUT_PFM) as handle:
        motif = motifs.read(handle, "jaspar")

        # Calculate pseudocounts
        motif.pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)  # type: ignore

        # Caclulate PWM
        pwm = list(map(list, zip(*[motif.pssm[nt] for nt in "ACGT"])))  # type: ignore

        # Write out PWM
        with open(OUTPUT_PWM, mode="w", encoding="utf-8") as op:
            for line in pwm:
                line = "\t".join(
                    ["{:7d}".format(round(j * 100)) for j in line]
                )  # type: ignore
                op.write(f"{line}\n")

# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()