from os import listdir
from snakemake.utils import min_version

# Parameters
ASSEMBLY = config["assembly"]
MATRICES = listdir(config["matrices"])
MATRIX_PROB = config["matrix_prob"]
MATRIX_SCAN = config["matrix_scan"]
CHROMOSOMES = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
MATRIX_DIR = config["matrices"]
# Settings
min_version("7.32.4")


rule all:
    input:
        expand(
            "results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}-sites.masked.genome.sorted.bed.gz",
            matrix=MATRICES,
            assembly=ASSEMBLY,
        ),
    default_target: True


rule compile_pwmscan:
    message:
        """
        Compiles PWMScan. Note c99 flag.
        Installs into resources/software by default.
        """
    input:
        prob=workflow.source_path(f"../../{MATRIX_PROB}"),
        scan=workflow.source_path(f"../../{MATRIX_SCAN}"),
    output:
        compiled_prob="resources/software/PWMScan/matrix_prob",
        compiled_scan="resources/software/PWMScan/matrix_scan",
    log:
        stdout="workflow/logs/compile_pwmscan.stdout",
        stderr="workflow/logs/compile_pwmscan.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        gcc -std=c99 -o {output.compiled_prob} {input.prob}
        gcc -std=c99 -o {output.compiled_scan} {input.scan}
        """

rule calculate_IntLogOdds:
    message:
        """
        Converts intput motif models to the IntLogOdds PWM format.
        """
    input:
        MATRIX_DIR + "/{matrix}",
    output:
        "results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}.IntLogOdds",
    conda:
        "../envs/tfbs-scan.yaml"
    log:
        stdout="workflow/logs/calculate_pwm_{matrix}_{assembly}.stdout",
        stderr="workflow/logs/calculate_pwm_{matrix}_{assembly}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    script:
        "../scripts/pwm/pwm.py"


rule calculate_probabilities:
    message:
        """
        Generates pval distribution for given PWM.
        """
    input:
        pwm=rules.calculate_IntLogOdds.output,
        matrix_prob=rules.compile_pwmscan.output.compiled_prob,
    output:
        temp("results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}-pvals.raw"),
    log:
        stdout="workflow/logs/calculate_probabilities_{matrix}_{assembly}.stdout",
        stderr="workflow/logs/calculate_probabilities_{matrix}_{assembly}.stderr",
    threads: 1
    conda:
        "../envs/tfbs-scan.yaml"
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
        "results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}-pvals.txt",
    log:
        stdout="workflow/logs/process_probabilities_{matrix}_{assembly}.stdout",
        stderr="workflow/logs/process_probabilities_{matrix}_{assembly}.stderr",
    threads: 1
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        """
        sed 's/%//g' {input} | sort -k1n > {output}
        """


rule calculate_cutoff:
    message:
        """
        Relative threshold is percentage of best PWM.
        For 80%, format as interger 80. Makes life easier.
        """
    input:
        rules.process_probabilities.output,
    output:
        "results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}-coeff.txt",
    params:
        pthresh=0.05,
        rthresh=80,
    log:
        stdout="workflow/logs/calculate_cutoff_{matrix}_{assembly}.stdout",
        stderr="workflow/logs/calculate_cutoff_{matrix}_{assembly}.stderr",
    threads: 1
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        """
        awk '{{if($2 < {params.pthresh} && $3>={params.rthresh}) print $1}}' {input} | head -n1 > {output}
        """


rule scan_chromosome:
    message:
        """
        Scans reference genome for matches againts input motif.
        Output compressed.
        """
    input:
        pwm=rules.calculate_IntLogOdds.output,
        cut=rules.calculate_cutoff.output,
        ref="results/tfbs-scan/{assembly}/mask/{chrom}.masked.fa",
        matrix_scan=rules.compile_pwmscan.output.compiled_scan,
    output:
        temp(
            "results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}-sites.masked.{chrom}.bed"
        ),
    log:
        stdout="workflow/logs/scan_chromosome_{matrix}_{assembly}_{chrom}.stdout",
        stderr="workflow/logs/scan_chromosome_{matrix}_{assembly}_{chrom}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        {input.matrix_scan} -m {input.pwm} -c $(cat {input.cut}) {input.ref} > {output}
        """


rule assemble_scan:
    message:
        """
        Assembles individual chromosome scans into genome sites file.
        """
    input:
        expand(
            "results/tfbs-scan/{assembly}/scan/{{matrix}}/{{matrix}}-sites.masked.{chrom}.bed",
            assembly=ASSEMBLY,
            chrom=CHROMOSOMES,
        ),
    output:
        "results/tfbs-scan/{assembly}/scan/{matrix}/{matrix}-sites.masked.genome.sorted.bed.gz",
    log:
        stdout="workflow/logs/assemble_scan_{assembly}_{matrix}_hg38.stdout",
        stderr="workflow/logs/assemble_scan_{assembly}_{matrix}_hg38.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        cat {input} | sort -k 1,1 -k2,2n | gzip > {output}
        """
