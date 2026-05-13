## What Changed vs HW2

### `main.nf`

- Removed hardcoded `conda` paths from every process — dependencies are now managed entirely via `nextflow.config`
- Replaced hardcoded `memory '8 GB'` with `label 'process_medium'` — resource requirements are now defined centrally per profile
- Added `tag` directive to all processes for readable logs
- Added `emit: bai` to `mapping_wf` — required to pass BAM index to variant calling
- Added `include { BCFTOOLS_MPILEUP }` — imports nf-core bcftools/mpileup module 
- Added `variant_calling_wf` — runs bcftools mpileup | call on the BAM file 
- Added `vcf` to `publish` and `output` blocks → results saved to `results/vcf/`

### `nextflow.config` — created from scratch

- All pipeline parameters (`input_reads_folder`, `sra_id`, `reference`, `threads`, `save_mpileup`)
- Three execution profiles
- Per-process container images for the `container` profile
- Automatic retry on OOM errors (exit codes 137, 143, 139)

## Configuration Profiles

| Profile     | Executor | Dependency management                            |
| ----------- | -------- | ------------------------------------------------ |
| `local`     | local    | Paths to pre-built conda envs on current machine |
| `cluster`   | cluster  | Per-process conda `.yml` files are created       |
| `container` | local    | Docker images pulled from Docker Hub             |

## nf-core Module

`bcftools/mpileup` was installed using nf-core tools:
```bash
nf-core modules install bcftools/mpileup
```
The module runs `bcftools mpileup | bcftools call` on the sorted BAM produced by `map_reads`, outputting a compressed VCF with index.


## Docker Images

Six images built manually and pushed to Docker Hub (`egoncharov/*`):

| Image | Tool |
|---|---|
| `egoncharov/sra-tools:3.0.10` | fasterq-dump |
| `egoncharov/fastqc:0.12.1` | FastQC |
| `egoncharov/trimmomatic:0.39` | Trimmomatic |
| `egoncharov/spades:3.15.5` | SPAdes assembler |
| `egoncharov/bowtie2-samtools:2.5.1` | Bowtie2 + Samtools |
| `egoncharov/matplotlib:3.8` | Python + Matplotlib |

![[Pasted image 20260513234831.png]]

`biocontainers/bcftools:1.18` (used by the nf-core module) is an official image — no custom build required.

## Usage

```bash
# Local (pre-built conda envs)
nextflow run main.nf -profile local --input_reads_folder /data/fastq

# Cluster (conda from .yml)
nextflow run main.nf -profile cluster --sra_id SRR12345678

# Container (Docker)
nextflow run main.nf -profile container --input_reads_folder /data/fastq

# Resume after failure
nextflow run main.nf -profile container --input_reads_folder /data/fastq -resume
```

If `--reference` is not provided, SPAdes assembles a reference de novo from trimmed reads. The resulting `scaffolds.fasta` is used for both read mapping and variant calling.