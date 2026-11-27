#!/usr/bin/env nextflow

/*
 * fastqgwas: a small genome-wide association study from raw reads. Aligns each sample's
 * reads, calls SNPs in all samples together, and tests each SNP for association with
 * case/control status. Run with --help to see the parameters.
 */

def helpText() {
    return """
    Usage:
      nextflow run ss-lab-imperial/fastqgwas -profile imperial \\
          --input samplesheet.csv --outdir nf_out -params-file gwas.yml

    Parameters:
      --input       samplesheet: CSV with columns sample,fastq_1, one FASTQ per line
                    (absolute paths, or relative to the samplesheet's folder)          (required)
      --fasta       reference genome, FASTA                                             (required)
      --phenotypes  #IID and PHENO columns; PHENO 2 = case, 1 = control                 (required)
      --outdir      folder for the final files                                          (required)
      --max_reads   reads aligned per sample (the first N in each file)  (default: ${params.max_reads})
      --min_qual    lowest variant QUAL kept                             (default: ${params.min_qual})

    Profiles:
      -profile imperial   Imperial's HPC: each step is a PBS job, tools come from containers
      -profile test       the small simulated dataset in test_data/ (add it to another
                          profile, e.g. -profile test,imperial; --outdir still needed)

    Outputs, in --outdir:
      calls.vcf.gz   biallelic SNPs called jointly in all samples
      gwas.tsv       one row per SNP: chrom, pos, id, ref, alt, n, odds ratio, p
      pipeline_info/ execution trace (.txt), report and timeline (.html) and step diagram
    """.stripIndent()
}

// Gives each sample its own copy of the reference, ref.fa, so that no two samples share files,
// and indexes it: for bowtie2 (ref.*.bt2) and for variant calling (ref.fa.fai)
process PREP_REF {
    label 'process_single'
    tag "${sample}"

    input:
    tuple val(sample), path(fasta, stageAs: 'input.fa')

    output:
    tuple val(sample), path('ref*'), emit: ref

    script:
    """
    cp input.fa ref.fa
    bowtie2-build --quiet --threads ${task.cpus} ref.fa ref
    samtools faidx ref.fa
    """
}

// Aligns one sample's reads to that sample's copy of the reference and sorts the alignments by
// position
process ALIGN {
    label 'process_single'
    tag "${sample}"

    input:
    tuple val(sample), path(fastq), path(reference)

    output:
    path "${sample}.bam", emit: bams
    path "${sample}.bam.bai", emit: bais

    script:
    """
    set -o pipefail
    bowtie2 -p ${task.cpus} -u ${params.max_reads} -x ref -U ${fastq} --rg-id ${sample} --rg SM:${sample} 2> ${sample}.log \\
        | samtools sort -o ${sample}.bam
    samtools index ${sample}.bam
    """
}

// Joint calling: every sample's alignments go into one mpileup, so all samples get a genotype
// at every site
process CALL {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path bams
    path bais
    path reference

    output:
    path 'calls.vcf.gz', emit: vcf

    script:
    """
    set -o pipefail
    ls *.bam > bams.txt
    bcftools mpileup --threads ${task.cpus} -f ref.fa -b bams.txt -a AD,DP -Ou \
        | bcftools call -mv -Ou \
        | bcftools view -m2 -M2 -v snps -i 'QUAL >= ${params.min_qual}' -Oz -o calls.vcf.gz
    """
}

// Logistic regression of case/control status on each SNP's genotype (plink2 --glm)
process ASSOC {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    path vcf
    path phenotypes

    output:
    path 'gwas.tsv'

    script:
    """
    plink2 --vcf ${vcf} --pheno ${phenotypes} --allow-extra-chr --set-all-var-ids '@:#' \
        --glm allow-no-covars --out gwas
    # Keep the useful columns (found by name), dropping SNPs plink2 could not test
    awk -F'\\t' 'BEGIN { OFS = "\\t" }
        NR == 1 { for (i = 1; i <= NF; i++) c[\$i] = i
                  print "chrom", "pos", "id", "ref", "alt", "n", "or", "p"; next }
        \$c["P"] != "NA" && (!("ERRCODE" in c) || \$c["ERRCODE"] == ".") {
            print \$c["#CHROM"], \$c["POS"], \$c["ID"], \$c["REF"], \$c["ALT"], \$c["OBS_CT"], \$c["OR"], \$c["P"] }' \
        gwas.PHENO.glm.logistic* > gwas.tsv
    """
}

workflow {
    if (params.help) {
        log.info helpText()
        return
    }
    if (!params.input || !params.fasta || !params.phenotypes || !params.outdir) {
        error "Set --input, --fasta, --phenotypes and --outdir (see --help)"
    }

    def sheet = file(params.input)
    if (!sheet.exists()) {
        error "No samplesheet at ${params.input}"
    }
    def header = sheet.readLines()[0].split(',')*.trim()
    if (!('sample' in header && 'fastq_1' in header)) {
        error "The samplesheet's first line must name the columns sample,fastq_1 (it is: ${header.join(',')})"
    }
    reads = channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            def fastq = sheet.parent.resolve(row.fastq_1)
            if (!fastq.exists()) {
                error "Sample ${row.sample}: no FASTQ file at ${fastq}. Relative paths in the samplesheet are read from the samplesheet's folder (${sheet.parent}), not from where you run nextflow"
            }
            tuple(row.sample, fastq)
        }

    // One copy of the reference per sample, then each sample's reads with its own copy
    PREP_REF(reads.map { sample, fastq -> tuple(sample, file(params.fasta)) })
    ALIGN(reads.join(PREP_REF.out.ref))
    // Variant calling needs only one copy: any sample's will do
    CALL(ALIGN.out.bams.collect(), ALIGN.out.bais.collect(), PREP_REF.out.ref.first().map { sample, reference -> reference })
    ASSOC(CALL.out.vcf, file(params.phenotypes))
}
