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
        extra=PERBASE_EXTRA
    threads: 24
    shell:
        """
        perbase base-depth {params.extra} --threads {threads} \
            --reference_fasta "{input.ref}" -Z "{input.bam}" -o "{output}" \
            1> "{log}" 2>&1
        """


""" 
TO-DO: Parse perbase output into a cleaner, annotated table that filters out unnecessary info,
summarizes the important stuff, and adds gene-level annotation information (chromsome, strand)
"""
