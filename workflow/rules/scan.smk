from snakemake.utils import min_version

# Parameters
ASSEMBLY = config["assembly"]
PROFILES = [i.split("|")[1] for i in config["targets"]]

# Settings
min_version("7.32.4")


# WC constraints - JASPAR matrix format
wildcard_constraints:
    PROFILE="[a-zA-Z\d]{6}.{1}\d{1}",


rule all:
    input:
        expand(
            "results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-sites.masked.bed.gz",
            ASSEMBLY=ASSEMBLY,
            PROFILE=PROFILES,
        ),


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
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
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
        "results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-pwm.txt",
    conda:
        "../envs/tfbs-scan.yaml"
    log:
        stdout="workflow/logs/calculate_pwm_{PROFILE}_{ASSEMBLY}.stdout",
        stderr="workflow/logs/calculate_pwm_{PROFILE}_{ASSEMBLY}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
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
        temp("results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-pvals.raw"),
    log:
        stdout="workflow/logs/calculate_probabilities_{PROFILE}_{ASSEMBLY}.stdout",
        stderr="workflow/logs/calculate_probabilities_{PROFILE}_{ASSEMBLY}.stderr",
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
        pvals="results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-pvals.txt",
        coeff="results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-coeff.txt",
    params:
        pthresh=0.05,
        rthresh=80,
    log:
        stdout="workflow/logs/process_probabilities_{PROFILE}_{ASSEMBLY}.stdout",
        stderr="workflow/logs/process_probabilities_{PROFILE}_{ASSEMBLY}.stderr",
    threads: 1
    conda:
        "../envs/tfbs-scan.yaml"
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
        "resources/data/genome/{ASSEMBLY}/{ASSEMBLY}.fa.gz",
    output:
        temp("results/tfbs-scan/{ASSEMBLY}/{ASSEMBLY}.fa"),
    log:
        stdout="workflow/logs/decompress_genome_{ASSEMBLY}.stdout",
        stderr="workflow/logs/decompress_genome_{ASSEMBLY}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        gunzip {input} -c > {output}
        """


rule mask_regions:
    input:
        genome=rules.decompress_genome.output,
        blacklist="resources/data/genome/{ASSEMBLY}/{ASSEMBLY}.blacklist.bed",
        exons="results/gencode/{ASSEMBLY}/gencode.{ASSEMBLY}.exons.protein_coding.bed",
    output:
        temp("results/tfbs-scan/{ASSEMBLY}/{ASSEMBLY}.masked_regions.bed"),
    params:
        gaps="resources/data/genome/{ASSEMBLY}/{ASSEMBLY}.gaps.bed",
    log:
        stdout="workflow/logs/mask_regions_{ASSEMBLY}.stdout",
        stderr="workflow/logs/mask_regions_{ASSEMBLY}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        cat {params.gaps} {input.exons} {input.blacklist} |
        vawk '{{print $1, $2, $3}}' | 
        vawk '!seen[$1, $2, $3]++' |
        sort -k 1,1 -k2,2n > {output}
        """


rule mask_genome:
    input:
        regions=rules.mask_regions.output,
        genome=rules.decompress_genome.output,
    output:
        masked_genome=temp("results/tfbs-scan/{ASSEMBLY}/hg38/hg38.custom-mask.fa"),
    log:
        stdout="workflow/logs/mask_genome_{ASSEMBLY}.stdout",
        stderr="workflow/logs/mask_genome_{ASSEMBLY}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        bedtools maskfasta -fi {input.genome} -bed {input.regions} -fo {output.masked_genome}
        """


rule scan_genome:
    message:
        """
        Scans reference genome for matches againts input motif.
        Output compressed.
        """
    input:
        pwm=rules.calculate_pwm.output,
        cut=rules.process_probabilities.output.coeff,
        ref=rules.mask_genome.output.masked_genome,
        matrix_scan="resources/software/PWMScan/matrix_scan",
    output:
        "results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-sites.masked.bed.gz",
    log:
        stdout="workflow/logs/scan_genome_{PROFILE}_{ASSEMBLY}.stdout",
        stderr="workflow/logs/scan_genome_{PROFILE}_{ASSEMBLY}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        {input.matrix_scan} -m {input.pwm} -c $(cat {input.cut}) {input.ref} |
        gzip > {output}
        """


rule bed2fasta:
    message:
        """
        Some utilities want FASTA format, useful to have on hand.
        - Note -s, keep strand information correct. Reports rev comp for neg seqs.
        Output compressed.
        """
    input:
        genome=rules.decompress_genome.output,
        sites=rules.scan_genome.output,
    output:
        "results/tfbs-scan/{ASSEMBLY}/{PROFILE}/{PROFILE}-sites.masked.fa.gz",
    log:
        stdout="workflow/logs/bed2fasta_{PROFILE}_{ASSEMBLY}.stdout",
        stderr="workflow/logs/bed2fasta_{PROFILE}_{ASSEMBLY}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        bedtools getfasta -fi {input.genome} -bed {input.sites} -name -s |
        gzip > {output}
        """
