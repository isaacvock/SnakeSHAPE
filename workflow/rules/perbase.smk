### Need to index BAM file first
rule index_transcriptome_bam:
    input:
        "results/align/{sample}.transcriptome.bam",
    output:
        "results/align/{sample}.transcriptome.bam.bai",
    conda:
        "../envs/genomictools.yaml"
    log:
        "logs/index_transcriptome_bam/{sample}.log",
    threads: 8
    shell:
        """
        samtools index -@ {threads} {input} 1> {log} 2>&1
        """

### Need to index FASTA file first
rule index_rsem_transcriptome_fasta:
    input:
        RSEM_TRANSCRIPTOME_FASTA,
    output:
        RSEM_TRANSCRIPTOME_FASTA_FAI,
    log:
        "logs/perbase/rsem_transcriptome_faidx.log",
    conda:
        "../envs/genomictools.yaml"
    threads: 2
    shell:
        """
        samtools faidx "{input}" 1> "{log}" 2>&1
        """


### Get per-nucleotide data
rule perbase:
    input:
        bam="results/align/{sample}.transcriptome.bam",
        ref=RSEM_TRANSCRIPTOME_FASTA,
        fai=RSEM_TRANSCRIPTOME_FASTA_FAI,
    output:
        "results/perbase/{sample}.tsv.gz",
    log:
        "logs/perbase/{sample}.log",
    conda:
        "../envs/perbase.yaml"
    params:
        extra=PERBASE_EXTRA,
    threads: 24
    shell:
        """
        perbase base-depth {params.extra} --threads {threads} \
            --ref-fasta "{input.ref}" -Z "{input.bam}" -o "{output}" \
            1> "{log}" 2>&1
        """
