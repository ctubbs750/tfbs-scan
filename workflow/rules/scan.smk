from sys import path
from pathlib import Path
from pandas import read_csv, merge

# Parameters
PROFILES = config["TFBS-SCAN"]["TARGETS"]
ASSEMBLY = config["GENOME"]["ASSEMBLY"]


# WC constraints - JASPAR matrix format
wildcard_constraints:
    PROFILE="[a-zA-Z\d]{6}.{1}\d{1}",


rule all:
    input:
        expand("results/tfbs-scan/{PROFILE}/{PROFILE}-sites.bed.gz", PROFILE=PROFILES),
        expand("results/tfbs-scan/{PROFILE}/{PROFILE}-sites.fa.gz", PROFILE=PROFILES),


rule download_jaspar:
    message:
        """
        Downloads motif matrix from the Jaspar database.
        """
    output:
        "resources/data/jaspar/{PROFILE}.jaspar",
    params:
        target="https://jaspar.elixir.no/api/v1/matrix/{PROFILE}/?format=jaspar",
    log:
        stdout="workflow/logs/download_jaspar_{PROFILE}.stdout",
        stderr="workflow/logs/download_jaspar_{PROFILE}.stderr",
    threads: 2
    shell:
        """
        curl {params.target} -o {output}
        """


rule calculate_pwm:
    """
    Converts counts matrix to weight matrix.
    """
    input:
        rules.download_jaspar.output,
    output:
        "results/tfbs-scan/{PROFILE}/{PROFILE}-pwm.txt",
    conda:
        "../envs/tfbs-scan.yaml"
    log:
        stdout="workflow/logs/calculate_pwm_{PROFILE}.stdout",
        stderr="workflow/logs/calculate_pwm_{PROFILE}.stderr",
    threads: 2
    script:
        "../scripts/pwm.py"


rule calculate_probabilities:
    message:
        """
        Generates pval distribution for given PWM.
        """
    input:
        pwm=rules.calculate_pwm.output,
        matrix_prob="resources/software/PWMScan/matrix_prob",
    output:
        temp("results/tfbs-scan/{PROFILE}/{PROFILE}-pvals.raw"),
    log:
        stdout="workflow/logs/calculate_probabilities_{PROFILE}.stdout",
        stderr="workflow/logs/calculate_probabilities_{PROFILE}.stderr",
    threads: 2
    shell:
        """
        {input.matrix_prob} {input} > {output}
        """


rule process_probabilities:
    message:
        """
        Relative threshold is percentage of best PWM.
        For 80%, format as interger 80. Makes life easier.
        """
    input:
        rules.calculate_probabilities.output,
    output:
        pvals="results/tfbs-scan/{PROFILE}/{PROFILE}-pvals.txt",
        coeff="results/tfbs-scan/{PROFILE}/{PROFILE}-coeff.txt",
    params:
        pthresh=0.05,
        rthresh=80,
    log:
        stdout="workflow/logs/process_probabilities_{PROFILE}.stdout",
        stderr="workflow/logs/process_probabilities_{PROFILE}.stderr",
    threads: 2
    shell:
        """
        set +o pipefail;
        sed 's/%//g' {input} | 
        sort -k1n > {output.pvals}
        awk '{{if($2 < {params.pthresh} && $3>={params.rthresh}) print $1}}' {output.pvals} |
        head -n1 > {output.coeff}
        """


rule decompress_genome:
    message:
        """
        Necessary to scan uncompressed genome.
        """
    input:
        f"resources/data/genome/{ASSEMBLY}/{ASSEMBLY}.fa.gz",
    output:
        f"resources/data/genome/{ASSEMBLY}/{ASSEMBLY}.fa",
    log:
        stdout=f"workflow/logs/decompress_genome_{ASSEMBLY}.stdout",
        stderr=f"workflow/logs/decompress_genome_{ASSEMBLY}.stderr",
    threads: 2
    shell:
        """
        gunzip {input}
        """


rule scan_genome:
    message:
        """
        Scans reference genome for matches againts input motif.
        """
    input:
        pwm=rules.calculate_pwm.output,
        cut=rules.process_probabilities.output.coeff,
        ref=rules.decompress_genome.output,
        matrix_scan="resources/software/PWMScan/matrix_scan",
    output:
        "results/tfbs-scan/{PROFILE}/{PROFILE}-sites.bed",
    log:
        stdout="workflow/logs/scan_genome_{PROFILE}.stdout",
        stderr="workflow/logs/scan_genome_{PROFILE}.stderr",
    threads: 4
    shell:
        """
        {input.matrix_scan} -m {input.pwm} -c $(cat {input.cut}) {input.ref} > {output}
        """


rule bed2fasta:
    message:
        """
        Some utilities want FASTA format, useful to have on hand.
        - Note -s, keep strand information correct. Reports rev comp for neg seqs.
        """
    input:
        genome=rules.decompress_genome.output,
        sites=rules.scan_genome.output,
    output:
        "results/tfbs-scan/{PROFILE}/{PROFILE}-sites.fa",
    conda:
        "../envs/tfbs-scan.yaml"
    log:
        stdout="workflow/logs/bed2fasta_{PROFILE}.stdout",
        stderr="workflow/logs/bed2fasta_{PROFILE}.stderr",
    threads: 4
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.sites} -name -s > {output}
        """


rule compress_sites:
    message:
        """
        Compress sites files to minimize storage.
        """
    input:
        bed=rules.scan_genome.output,
        fa=rules.bed2fasta.output,
    output:
        "results/tfbs-scan/{PROFILE}/{PROFILE}-sites.bed.gz",
        "results/tfbs-scan/{PROFILE}/{PROFILE}-sites.fa.gz",
    log:
        stdout="workflow/logs/compress_sites_{PROFILE}.stdout",
        stderr="workflow/logs/compress_sites_{PROFILE}.stderr",
    threads: 4
    shell:
        """
        gzip {input.bed} && gzip {input.fa}
        """
