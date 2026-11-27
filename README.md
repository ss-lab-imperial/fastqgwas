# fastqgwas

A small genome-wide association study (GWAS) from raw reads: align each sample's reads to the
reference, call SNPs in all samples together, and test each SNP for association with
case/control status.

## Quick start

On Imperial's HPC, 

```bash
module load Nextflow
nextflow run ss-lab-imperial/fastqgwas -r 1.0.0 -profile imperial \
    --outdir nf_out \
    --input samplesheet.csv \
    -params-file gwas.yml \
    -work-dir /rds/general/user/${USER}/ephemeral/nf_work_gwas \
    -resume
```

with `gwas.yml` holding the other inputs:

```yaml
fasta: <path_to_your_data>/ref.fa
phenotypes: <path_to_your_data>/phenotypes.tsv
```

To try it on the small test dataset in `test_data/` first:

```bash
nextflow run ss-lab-imperial/fastqgwas -r 1.0.0 -profile test,imperial --outdir test_out
```

(`-profile test` supplies `--input` and the `-params-file` inputs; `--outdir` is always needed.)



## Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `input` | CSV with columns `sample,fastq_1` (absolute paths, or relative to the CSV) | required |
| `fasta` | Reference genome (FASTA) | required |
| `phenotypes` | Tab-separated, columns `#IID` and `PHENO` (2 = case, 1 = control) | required |
| `outdir` | Folder for the final files | required |
| `max_reads` | Reads aligned per sample: the first N in each FASTQ file (the reads are in random order, so these are a random sample). Your FASTQ files are larger than a GWAS needs. | `300000` |
| `min_qual` | Lowest variant QUAL kept | `20` |

## Processes

| Process | What it does | Tools |
| --- | --- | --- |
| `PREP_REF` | Prepares the reference | samtools |
| `ALIGN` | Aligns the reads to the reference and sorts the alignments | bowtie2, samtools |
| `CALL` | Calls SNPs in all samples together | bcftools |
| `ASSOC` | Logistic regression of case/control status on each SNP | plink2 |

bowtie2 needs an index of the reference, which the workflow builds with `bowtie2-build`.

Tools come from containers (the images nf-core's own modules use), which `-profile imperial` runs with Singularity.

## Outputs

- `calls.vcf.gz`: biallelic SNPs with a genotype for every sample
- `gwas.tsv`: one row per SNP: `chrom pos id ref alt n or p`
- `pipeline_info/`: how the run went, one set of files per run (named by date and time):
  - `execution_trace_*.txt`: one row per task, with its status, exit code, run time, CPU and memory
  - `execution_report_*.html`: summary of the run, with resource use per step
  - `execution_timeline_*.html`: when each task started and finished
  - `pipeline_dag_*.html`: diagram of the steps
