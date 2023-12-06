"""Main ruleset for the pwm-scan module"""

from sys import path
from pathlib import Path
from pandas import read_csv, merge

# I/O
INSTALL_DIR = config["INSTALL_DIR"]
PROCESS_DIR = config["PROCESS_DIR"]

# Parameters
PROFILES = config["TARGETS"]

# Software
MATRIX_PROB = config["PWMSCAN"]["MATRIX_PROB"]
MATRIX_SCAN = config["PWMSCAN"]["MATRIX_SCAN"]


# WC constraints - JASPAR matrix format
wildcard_constraints:
    PROFILE="[a-zA-Z\d]{6}.{1}\d{1}",


rule all:
    input:
        expand(PROCESS_DIR + "/{PROFILE}/{PROFILE}-sites.bed.gz", PROFILE=PROFILES),
        expand(PROCESS_DIR + "/{PROFILE}/{PROFILE}-sites.fa.gz", PROFILE=PROFILES),


rule compile_pwmscan:
    message:
        """
        c99 mode
        """
    input:
        prob=workflow.source_path(f"{MATRIX_PROB}.c"),
        scan=workflow.source_path(f"{MATRIX_SCAN}.c"),
    output:
        compiled_prob=MATRIX_PROB,
        compiled_scan=MATRIX_SCAN,
    shell:
        """
        gcc -std=c99 -o matrix_prob {input.prob} &&
        gcc -std=c99 -o matrix_scan {input.scan}
        """


rule download_jaspar:
    message:
        """
        - Downloads motif PCM from the Jaspar database.
        """
    output:
        INSTALL_DIR + "/{PROFILE}.jaspar",
    params:
        target="https://jaspar.elixir.no/api/v1/matrix/{PROFILE}/?format=jaspar",
    shell:
        "curl {params.target} -o {output}"


rule calculate_pwm:
    """
    - Converts counts matrix to weight matrix.
    """
    input:
        rules.download_jaspar.output,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-pwm.txt",
    conda:
        "genome"
    script:
        "../scripts/pwm.py"


rule calc_probabilities:
    message:
        """
        - Generates pval distribution for given PWM.
        """
    input:
        rules.calculate_pwm.output,
    output:
        temp(PROCESS_DIR + "/{PROFILE}/{PROFILE}-pvals.raw"),
    params:
        matrix_prob=rules.compile_pwmscan.output.compiled_prob,
    shell:
        """
        {params.matrix_prob} {input} > {output}
        """


rule process_probabilities:
    message:
        """
        - Relative threshold is percentage of best PWM.
        For 80%, format as interger 80. Makes life easier.
        """
    input:
        rules.calc_probabilities.output,
    output:
        pvals=PROCESS_DIR + "/{PROFILE}/{PROFILE}-pvals.txt",
        coeff=PROCESS_DIR + "/{PROFILE}/{PROFILE}-coeff.txt",
    params:
        pthresh=0.05,
        rthresh=80,
    shell:
        """
        set +o pipefail
        sed 's/%//g' {input} > {output.pvals}
        awk '{{if($2 < {params.pthresh} && $3=={params.rthresh}) print $1}}' {output.pvals} > {output.coeff}
        """


rule decompress_genome:
    message:
        """
        - d
        """
    input:
        "../resources/data/genome/hg38/hg38.fa.gz",
    output:
        "../resources/data/genome/hg38/hg38.fa",
    shell:
        """
        gunzip {input}
        """


rule scan_genome:
    message:
        """
        - Scans rerference genome for matches againts input motif.
        """
    input:
        matrix=rules.calculate_pwm.output,
        cutoff=rules.process_probabilities.output.coeff,
        genome=rules.decompress_genome.output,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-sites.bed",
    params:
        matrix_scan=rules.compile_pwmscan.output.compiled_scan,
    shell:
        """
        {params.matrix_scan} -m {input.matrix} -c $(cat {input.cutoff}) {input.genome} > {output}
        """


rule bed2fasta:
    message:
        """
        - Note -s, keep strand information correct. Reports rev comp for neg seqs.
        """
    conda:
        "../envs/pwmscan.yaml"
    input:
        genome=rules.decompress_genome.output,
        sites=rules.scan_genome.output,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-sites.fa",
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.sites} -name -s > {output}
        """


rule compress_sites:
    message:
        """
        - d
        """
    input:
        bed=rules.scan_genome.output,
        fa=rules.bed2fasta.output,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-sites.bed.gz",
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-sites.fa.gz",
    shell:
        """
        gzip {input.bed} && gzip {input.fa}
        """
