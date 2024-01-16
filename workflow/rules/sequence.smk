from snakemake.utils import min_version

# Parameters
ASSEMBLY = config["assembly"]
CHROMOSOMES = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"]

# Settings
min_version("7.32.4")


rule all:
    input:
        expand(
            "results/tfbs-scan/{assembly}/mask/{chrom}.masked.fa",
            assembly=ASSEMBLY,
            chrom=CHROMOSOMES,
        ),
    default_target: True


rule decompress_genome:
    message:
        """
        Necessary to scan genome.
        """
    input:
        "resources/data/genome/{assembly}/{assembly}.fa.gz",
    output:
        temp("results/tfbs-scan/{assembly}/{assembly}.fa"),
    log:
        stdout="workflow/logs/decompress_genome_{assembly}.stdout",
        stderr="workflow/logs/decompress_genome_{assembly}.stderr",
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
        blist="resources/data/genome/{assembly}/{assembly}.blacklist.bed",
        exons="results/gencode/{assembly}/gencode.{assembly}.exons.protein_coding.bed",
    output:
        temp("results/tfbs-scan/{assembly}/mask/{assembly}.masked_regions.bed"),
    params:
        gaps="resources/data/genome/{assembly}/{assembly}.gaps.bed",
    log:
        stdout="workflow/logs/mask_regions_{assembly}.stdout",
        stderr="workflow/logs/mask_regions_{assembly}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    threads: 1
    shell:
        """
        cat {params.gaps} {input.exons} {input.blist} |
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
        temp("results/tfbs-scan/{assembly}/mask/hg38.custom-mask.fa"),
    log:
        stdout="workflow/logs/mask_genome_{assembly}.stdout",
        stderr="workflow/logs/mask_genome_{assembly}.stderr",
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
        """
    input:
        rules.mask_genome.output,
    output:
        "results/tfbs-scan/{assembly}/mask/{chrom}.masked.fa",
    params:
        chromosome=lambda wc: wc.chrom,
    log:
        stdout="workflow/logs/split_genome_{assembly}_{chrom}.stdout",
        stderr="workflow/logs/split_genome_{assembly}_{chrom}.stderr",
    conda:
        "../envs/tfbs-scan.yaml"
    # High thread-count seems to be necessary to avoid issues with parallel access
    threads: 3
    shell:
        """
        faidx --regex "{params.chromosome}" --out {output} {input}
        """
