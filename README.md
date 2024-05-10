# DataForAsobi

## はじめに

- ここはPlayground(遊び場）である
- 現実世界のデータはサイズが大きすぎるため、繰り返しの試行に向いていない。このことが学習を困難にしている。
- そこで小さな仮想のデータを作って、シミュレーションを行うことで勉強することを試みる。
- 参考資料 [Biostar Handbook](https://www.biostarhandbook.com/)

## リファレンスゲノム

まずは小さな架空のリファレンスゲノムを作ってみる。

```python
import random

# 染色体名と長さを定義
chromosomes = {
    'chr1': 2000,
    'chr2': 1800,
    'chr3': 1600,
    'chr4': 1400,
    'chr5': 1200,
}

# 塩基をランダムに生成する関数
def generate_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

# FASTAファイルを出力
for chromosome, length in chromosomes.items():
    sequence = generate_sequence(length)
    print(f">{chromosome}\n{sequence}")
```

```
python generate_reference.py > ref.fa
```

bgzip による圧縮

```
bgzip -k ref.fa
```

faidx を作る

```
samtools faidx ref.fa
```

bwa のインデックスを作る

```
bwa index ref.fa
```

2bit ファイルを作る

```
faToTwoBit ref.fa ref.2bit
```


## 突然変異とショートリードのシミュレーション

[wgsim.cr](https://github.com/kojix2/wgsim.cr) を使います

```
mkdir s1 s2 s3
```

mutationの発生

```
wgsim.cr mut -S 201 -i 0.0002 -d 0.0002 -I 0.5 -D 0.5 ref/ref.fa > s1/s1.fa 2> s1/s1.mutation.txt
wgsim.cr mut -S 301 -i 0.0002 -d 0.0002 -I 0.5 -D 0.5 ref/ref.fa > s2/s2.fa 2> s2/s2.mutation.txt
wgsim.cr mut -S 401 -i 0.0002 -d 0.0002 -I 0.5 -D 0.5 ref/ref.fa > s3/s3.fa 2> s3/s3.mutation.txt
```

シーケンスのシミュレーション

```
wgsim.cr seq -S 201 s1/s1.fa s1/s1_1.fa s1/s1_2.fa
wgsim.cr seq -S 301 s1/s1.fa s2/s2_1.fa s2/s2_2.fa
wgsim.cr seq -S 401 s1/s1.fa s3/s3_1.fa s3/s3_2.fa
```

## アラインメント

```snakemake
SAMPLES = ["S1", "S2", "S3"]
REFERENCE = "ref/ref.fa"

rule all:
    input:
        expand("{sample}/{sample}.sorted.bam", sample=[s.lower() for s in SAMPLES])

rule bwa_mem:
    input:
        ref=REFERENCE,
        r1="{sample}/{sample}_1.fa.gz",
        r2="{sample}/{sample}_2.fa.gz"
    output:
        bam="{sample}/{sample}.bam"
    params:
        rg="@RG\\tID:{sample}\\tSM:{sample}\\tLB:{sample}"
    shell:
        "bwa mem -R {params.rg} {input.ref} {input.r1} {input.r2} > {output.bam}"

rule samtools_sort:
    input:
        "{sample}/{sample}.bam"
    output:
        "{sample}/{sample}.sorted.bam"
    shell:
        "samtools sort -o {output} {input}"
```

```
snakemake --cores 3
```

みたいな感じ

## bcftools

bcftools mpileup -O v -f ref.fa s*.sorted.bam > s123_genotypes.vcf
bcftools call -vm -Ov s123_genotypes.vcf | bcftools norm -Ov -f ref.fa -d all - > s123_variants.vcf

## freebayes

## mosdepth

quantile すごい

```
mosdepth -q 0:1:4:100:200: s1 s1.sorted.bam
```
