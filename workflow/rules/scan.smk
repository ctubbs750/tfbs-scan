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
    "sites.masked.genome.sorted.bed",
)

# Logo plot
PROFILE_LOGO = os.path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "scan",
    "{tf_name}",
    "{profile}",
    "{dataset}",
    "logo.png",
)

# PWM IC
PROFILE_IC = os.path.join(
    OUTPUT_DIR,
    ASSEMBLY,
    "scan",
    "{tf_name}",
    "{profile}",
    "{dataset}",
    "IntLogOdds.IC.tsv",
)

# ------------- #
# Params        #
# ------------- #

wildcard_constraints:
    tf_name="\w+",
    profile="MA\d{4}\.\d",
    dataset="[^/]+",
    chrom="chr[0-9XYM]+"

TF_NAMES, PROFILES, DATASETS = glob_wildcards(
    os.path.join(PROFILES_DIR, "{tf_name}", "{profile}", "{dataset}" + f".{FORMAT}")
)

SCAN_THRESHOLDS = "results/unibind/biosample_thresholds.txt"
DICT_THRESHOLDS = dict(
    (line.split()[0], (float(line.split()[-1])))
    for line in open(SCAN_THRESHOLDS)
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
        expand(
            PROFILE_LOGO,
            zip,
            tf_name=TF_NAMES,
            profile=PROFILES,
            dataset=DATASETS,
        ),


rule decompress_genome:
    message:
        "Necessary to scan genome."
    input:
        GENOME_FILE,
    output:
        temp(f"{path.splitext(GENOME_FILE)[0]}.fa"),
    log:
        stdout="workflow/logs/decompress_genome.stdout",
        stderr="workflow/logs/decompress_genome.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        "gunzip {input} -c > {output}"


rule mask_regions:
    message:
        "Makes exclude regions in BED format for genome masking."
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
    shell:
        """
        # Concatenate gaps, exons, and blacklist files
        cat {params.gaps} {input.exons} {input.blacklist} |
        # Print the first three columns of each line
        vawk '{{print $1, $2, $3}}' | 
        # Remove duplicate lines
        vawk '!seen[$1, $2, $3]++' |
        # Sort the output by the first two columns
        sort -k 1,1 -k2,2n > {output}
        """


rule mask_genome:
    message:
        "Creates custom-masked genome for scanning against."
    input:
        regions=rules.mask_regions.output,
        genome=rules.decompress_genome.output,
    output:
        temp(GENOME_MASK),
    log:
        stdout="workflow/logs/mask_genome.stdout",
        stderr="workflow/logs/mask_genome.stderr",  # Corrected the log file name
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        """
        # Use bedtools to mask the genome with the regions
        bedtools maskfasta -fi {input.genome} -bed {input.regions} -fo {output}
        """


rule split_genome:
    message:
        "Splits masked genome by chromosome for parallel scanning."
    input:
        rules.mask_genome.output,
    output:
        CHROMOSOME_MASK,
    params:
        chrom=lambda wc: wc.chrom,
    log:
        stdout="workflow/logs/split_genome_{chrom}.stdout",
        stderr="workflow/logs/split_genome_{chrom}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 12  # High thread-count to avoid issues with parallel access
    shell:
        """
        # Use faidx to split the genome by chromosome
        faidx --regex "{params.chrom}" --out {output} {input}
        """


rule compile_pwmscan:
    message:
        "Compiles PWMScan. Note c99 flag. Installs into resources/software by default."
    input:
        # prob=workflow.source_path(f"../../{MATRIX_PROB}"),
        # scan=workflow.source_path(f"../../{MATRIX_SCAN}"),
        prob="../scripts/scan/matrix_prob.c",
        scan="../scripts/scan/matrix_scan.c",
    output:
        compiled_prob=COMPILED_PROB,
        compiled_scan=COMPILED_SCAN,
    log:
        stdout="workflow/logs/compile_pwmscan.stdout",
        stderr="workflow/logs/compile_pwmscan.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        """
        # Compile the prob and scan files with gcc
        gcc -std=c99 -o {output.compiled_prob} {input.prob}
        gcc -std=c99 -o {output.compiled_scan} {input.scan}
        """


rule calculate_IntLogOdds:
    message:
        "Converts input motif models to the IntLogOdds PWM format."
    input:
        PROFILE,
    output:
        PROFILE_IntLogOdds,
    log:
        stdout="workflow/logs/calculate_pwm_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/calculate_pwm_{tf_name}_{profile}_{dataset}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
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
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        "{input.matrix_prob} {input.pwm} > {output}"


rule process_probabilities:
    message:
        "Relative threshold is percentage of best PWM. For 80%, format as integer 80. Makes life easier."
    input:
        rules.calculate_probabilities.output,
    output:
        PVALS_CLEANED,
    log:
        stdout="workflow/logs/process_probabilities_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/process_probabilities_{tf_name}_{profile}_{dataset}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    shell:
        """
        # Remove percentage signs and sort the input
        sed 's/%//g' {input} | sort -k1n > {output}
        """


# rule calculate_cutoff:
#     message:
#         """
#         Relative threshold is percentage of best PWM.
#         For 80%, format as integer 80. Makes life easier.
#         Second awk grabs the head -n1, avoids pipefail
#         """
#     input:
#         rules.process_probabilities.output,
#     output:
#         CUTOFF,
#     params:
#         pthresh=0.05,
#         rthresh=80,
#     log:
#         stdout="workflow/logs/calculate_cutoff_{tf_name}_{profile}_{dataset}.stdout",
#         stderr="workflow/logs/calculate_cutoff_{tf_name}_{profile}_{dataset}.stderr",
#     conda:
#         "../envs/tfbs-scan.yaml"
#     shell:
#         """
#         # Use awk to filter the input and print the first line that matches the conditions
#         awk '{{if($2 < {params.pthresh} && $3>={params.rthresh}) print $1}}' {input} |
#         awk 'FNR == 1' > {output}
#         """


rule calculate_cutoff:
    message:
        """
        Get interger cutoff form unibind map
        """
    input:
        rules.calculate_probabilities.output,
    output:
        CUTOFF,
    params:
        biosample=lambda wc: wc.dataset,
    log:
        stdout="workflow/logs/calculate_cutoff_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/calculate_cutoff_{tf_name}_{profile}_{dataset}.stderr",
    run:
        # Get the threshold from the unibind map
        threshold = DICT_THRESHOLDS[wildcards.dataset]
        # Get score from prob that is closest to threshold
        with open(input[0], "r") as f:
            for line in f:
                score = float(line.split()[0])
                pval = float(line.split()[1])
                perc = float(line.split()[2])
                if perc < threshold:
                    cutoff = score
                    break
        # Write out
        with open(output[0], "w") as f:
            f.write(str(cutoff))

rule scan_chromosome:
    message:
        "Scans reference genome for matches against input motif. Output compressed."
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
    shell:
        """
        # Use matrix_scan to scan the reference genome for matches against the input motif
        {input.matrix_scan} -m {input.pwm} -c $(cat {input.cut}) {input.ref} > {output}
        """


rule assemble_scan:
    message:
        """
        Assembles individual chromosome scans into genome sites file.
        Vawk removes duplicate seqs on opposite strands due to palindromic seqs
        """
    input:
        expand(
            f"{OUTPUT_DIR}/{ASSEMBLY}/scan/"
            + "{{tf_name}}/{{profile}}/{{dataset}}/sites.masked.{chrom}.bed",
            tf_name=TF_NAMES, profile=PROFILES, chrom=CHROMOSOMES,
        ),
    output:
        temp(ASSEMBLED_SCAN),
    log:
        stdout="workflow/logs/assemble_scan_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/assemble_scan_{tf_name}_{profile}_{dataset}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 2
    shell:
        """
        cat {input} | vawk '!a[$1, $2, $3]++' | sort --parallel=2 -k 1,1 -k2,2n > {output}
        """


rule make_logo:
    message:
        "Saves PWM as logo and saves IC"
    input:
        PROFILE_IntLogOdds,
    output:
        logo=PROFILE_LOGO,
        ic=PROFILE_IC,
    params:
        dataset=lambda wc: wc.dataset,
    log:
        stdout="workflow/logs/make_logo_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/make_logo_{tf_name}_{profile}_{dataset}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    script:
        "../scripts/logo/logo.py"
