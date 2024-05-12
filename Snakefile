# Snakefile

import yaml

# Load configuration
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

# Load sample names
NORMAL_SAMPLES = config["samples"]["normal"]
TUMOR_SAMPLES = config["samples"]["tumor"]
REFERENCE = config["reference"]
GATK_IMAGE = config["gatk_image"]
DOCKER_DATA = config["docker_data"]
PAIRS = config["pairs"]

# Function to get relative path for Docker
def to_docker_path(path):
    return f"{DOCKER_DATA}/{path}"

# Default target
rule all:
    input:
        expand("{sample}/{sample}.sorted.bam.bai", sample=[s.lower() for s in NORMAL_SAMPLES + TUMOR_SAMPLES]),
        expand("mutect2/{normal}_{tumor}.vcf", zip, normal=[p['normal'].lower() for p in PAIRS], tumor=[p['tumor'].lower() for p in PAIRS])

rule bwa_mem:
    input:
        ref=REFERENCE,
        r1="{sample}/{sample}_1.fa",
        r2="{sample}/{sample}_2.fa"
    output:
        bam="{sample}/{sample}.bam"
    params:
        rg=lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}"
    shell:
        "bwa mem -R '{params.rg}' {input.ref} {input.r1} {input.r2} > {output.bam}"

rule samtools_sort:
    input:
        "{sample}/{sample}.bam"
    output:
        "{sample}/{sample}.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"

rule samtools_index:
    input:
        "{sample}/{sample}.sorted.bam"
    output:
        "{sample}/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule mutect2:
    input:
        normal_bam="{normal}/{normal}.sorted.bam",
        normal_bai="{normal}/{normal}.sorted.bam.bai",
        tumor_bam="{tumor}/{tumor}.sorted.bam",
        tumor_bai="{tumor}/{tumor}.sorted.bam.bai",
        ref=REFERENCE
    output:
        "mutect2/{normal}_{tumor}.vcf"
    run:
        docker_normal_bam = to_docker_path(input.normal_bam)
        docker_tumor_bam = to_docker_path(input.tumor_bam)
        docker_ref = to_docker_path(input.ref)
        docker_output = to_docker_path(output[0])

        shell(
            f"""
            docker run --rm -v $(pwd):{DOCKER_DATA} -u $(id -u):$(id -g) -it {GATK_IMAGE} \
                gatk Mutect2 \
                -R {docker_ref} \
                -I {docker_normal_bam} \
                -I {docker_tumor_bam} \
                --normal {wildcards.normal} \
                -O {docker_output}
            """
        )
