# Snakefile

import yaml

# Load configuration
with open("config.yml", "r") as f:
    config = yaml.safe_load(f)

# Load sample names
SAMPLES = config["samples"]
REFERENCE = config["reference"]

# Default target
rule all:
    input:
        expand("{sample}/{sample}.sorted.bam", sample=[s.lower() for s in SAMPLES])

rule bwa_mem:
    input:
        ref=REFERENCE,
        r1="{sample}/{sample}_1.fa",
        r2="{sample}/{sample}_2.fa"
    output:
        bam="{sample}/{sample}.bam"
    params:
        rg="@RG\\tID:{sample.upper()}\\tSM:{sample.upper()}\\tLB:{sample.upper()}"
    shell:
        "bwa mem -R '{params.rg}' {input.ref} {input.r1} {input.r2} > {output.bam}"

rule samtools_sort:
    input:
        "{sample}/{sample}.bam"
    output:
        "{sample}/{sample}.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"
