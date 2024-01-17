from os import listdir, path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

FORMAT = config["format"]
ASSEMBLY = config["assembly"]
OUTPUT_DIR = config["output_dir"]
GENOME_DIR = config["genome_dir"]
GENCODE_DIR = config["gencode_dir"]
MATRIX_PROB = config["matrix_prob"]
MATRIX_SCAN = config["matrix_scan"]
CHROMOSOMES = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]
PROFILES_DIR = config["profiles_dir"]

# ------------- #
# I/O           #
# ------------- #

# Input motif profiles
PROFILE = os.path.join(
    PROFILES_DIR, "{tf_name}", "{profile}", "{dataset}" + f".{FORMAT}"
)

# Motif profies IntLogOdds format
PROFILE_IntLogOdds = os.path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "scan",
    "{tf_name}",
    "{profile}",
    "{dataset}",
    "IntLogOdds.pwm",
)

# Genome fasta file
GENOME_FILE = os.path.join(GENOME_DIR, ASSEMBLY, f"{ASSEMBLY}.fa.gz")

# Encode blacklist
BLACKLIST = os.path.join(GENOME_DIR, ASSEMBLY, f"{ASSEMBLY}.blacklist.bed")

# Genome gaps
UCSC_GAPS = os.path.join(GENOME_DIR, ASSEMBLY, f"{ASSEMBLY}.gaps.bed")

# Gencode protein-coding exons
EXONS = os.path.join(GENCODE_DIR, ASSEMBLY, "gencode.hg38.exons.protein_coding.bed")

# Gencode protein-coding exons
MASK_REGIONS = os.path.join(
    OUTPUT_DIR, ASSEMBLY, "mask", f"{ASSEMBLY}.masked_regions.bed"
)

# Masked genome file
GENOME_MASK = os.path.join(OUTPUT_DIR, ASSEMBLY, "mask", f"{ASSEMBLY}.custom-mask.fa")

# Masked genome file split into chroms
CHROMOSOME_MASK = os.path.join(
    OUTPUT_DIR, ASSEMBLY, "mask", f"{ASSEMBLY}" + ".{chrom}.custom-mask.fa"
)

# Compiled matrix prob
COMPILED_PROB = os.path.join("resources", "software", "PWMScan", "matrix_prob")

# Compiled matrix scan
COMPILED_SCAN = os.path.join("resources", "software", "PWMScan", "matrix_scan")

# Raw pvals
PVALS_RAW = os.path.join(
    OUTPUT_DIR, ASSEMBLY, "scan", "{tf_name}", "{profile}", "{dataset}", "pvals_raw.txt"
)

# Clean pvals
PVALS_CLEANED = os.path.join(
    OUTPUT_DIR, ASSEMBLY, "scan", "{tf_name}", "{profile}", "{dataset}", "pvals.txt"
)

# Clean pvals
CUTOFF = os.path.join(
    OUTPUT_DIR, ASSEMBLY, "scan", "{tf_name}", "{profile}", "{dataset}", "cutoff.txt"
)

# Scanned chromosome
SCANNED_CHROMOSOME = os.path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "scan",
    "{tf_name}",
    "{profile}",
    "{dataset}",
    "sites.masked.{chrom}.bed",
)

# Scanned chromosome
ASSEMBLED_SCAN = os.path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "scan",
    "{tf_name}",
    "{profile}",
    "{dataset}",
    "sites.masked.genome.sorted.bed.gz",
)

# ------------- #
# Params        #
# ------------- #

TF_NAMES, PROFILES, DATASETS = glob_wildcards(
    os.path.join(PROFILES_DIR, "{tf_name}", "{profile}", "{dataset}" + f".{FORMAT}")
)

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        expand(
            ASSEMBLED_SCAN,
            zip,
            tf_name=TF_NAMES,
            profile=PROFILES,
            dataset=DATASETS,
        ),
    default_target: True


rule decompress_genome:
    message:
        """
        Necessary to scan genome.
        """
    input:
        GENOME_FILE,
    output:
        temp(f"{GENOME_FILE[:-3]}.fa"),
    log:
        stdout="workflow/logs/decompress_genome.stdout",
        stderr="workflow/logs/decompress_genome.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        gunzip {input} -c > {output}
        """


rule mask_regions:
    message:
        """
        Makes exclude regions in BED format for genome masking.
        """
    input:
        exons=EXONS,
        blacklist=BLACKLIST,
    output:
        temp(MASK_REGIONS),
    params:
        gaps=UCSC_GAPS,
    log:
        stdout="workflow/logs/mask_regions.stdout",
        stderr="workflow/logs/mask_regions.stderr",
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
    message:
        """
        Creates custom-masked genome for scanning against.
        """
    input:
        regions=rules.mask_regions.output,
        genome=rules.decompress_genome.output,
    output:
        temp(GENOME_MASK),
    log:
        stdout="workflow/logs/mask_genome.stdout",
        stderr="workflow/logs/mask_genomestderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        bedtools maskfasta -fi {input.genome} -bed {input.regions} -fo {output}
        """


rule split_genome:
    message:
        """
        Splits masked genome by chromsome from parallel scanning.
        High thread-count seems to be necessary to avoid issues with parallel access.
        Needs to be sufficiently high to scale with high number of cores.
        """
    input:
        rules.mask_genome.output,
    output:
        CHROMOSOME_MASK,
    params:
        chromosome=lambda wc: wc.chrom,
    log:
        stdout="workflow/logs/split_genome_{chrom}.stdout",
        stderr="workflow/logs/split_genome_{chrom}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 12
    shell:
        """
        faidx --regex "{params.chromosome}" --out {output} {input}
        """


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
        compiled_prob=COMPILED_PROB,
        compiled_scan=COMPILED_SCAN,
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
        PROFILE,
    output:
        PROFILE_IntLogOdds,
    conda:
        "../envs/tfbs-scan.yaml"
    log:
        stdout="workflow/logs/calculate_pwm_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/calculate_pwm_{tf_name}_{profile}_{dataset}.stderr",
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
        temp(PVALS_RAW),
    log:
        stdout="workflow/logs/calculate_probabilities_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/calculate_probabilities_{tf_name}_{profile}_{dataset}.stderr",
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
        PVALS_CLEANED,
    log:
        stdout="workflow/logs/process_probabilities_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/process_probabilities_{tf_name}_{profile}_{dataset}.stderr",
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
        CUTOFF,
    params:
        pthresh=0.05,
        rthresh=80,
    log:
        stdout="workflow/logs/calculate_cutoff_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/calculate_cutoff_{tf_name}_{profile}_{dataset}.stderr",
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
        ref=rules.split_genome.output,
        matrix_scan=rules.compile_pwmscan.output.compiled_scan,
    output:
        temp(SCANNED_CHROMOSOME),
    log:
        stdout="workflow/logs/scan_chromosome_{tf_name}_{profile}_{dataset}_{chrom}.stdout",
        stderr="workflow/logs/scan_chromosome_{tf_name}_{profile}_{dataset}_{chrom}.stderr",
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
            f"{OUTPUT_DIR}/{ASSEMBLY}/scan/"
            + "{{tf_name}}/{{profile}}/{{dataset}}/sites.masked.{chrom}.bed",
            chrom=CHROMOSOMES,
        ),
    output:
        ASSEMBLED_SCAN,
    log:
        stdout="workflow/logs/assemble_scan_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/assemble_scan_{tf_name}_{profile}_{dataset}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 2
    shell:
        """
        cat {input} | sort --parallel=2 -k 1,1 -k2,2n | gzip > {output}
        """