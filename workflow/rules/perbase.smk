### Get per-nucleotide data
rule perbase:
    input:
        bam="results/align/{sample}.transcriptome.bam",
        ref=RSEM_TRANSCRIPTOME_FASTA,
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
