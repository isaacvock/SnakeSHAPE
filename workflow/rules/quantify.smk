### Quantify gene expression (RNA-seq)
rule quantify_genes:
    input:
        bam="results/filter/{sample}.bam",
        gtf=config.get("annotation"),
    output:
        "results/quantify_genes/{sample}_gene.featureCounts",
    conda:
        "../envs/genomictools.yaml"
    params:
        strand=FEATURECOUNTS_STRAND,
        paired_end=lambda wc: "-p" if is_paired_end(wc.sample) else "",
        extra=FEATURECOUNTS_EXTRA,
    log:
        "logs/quantify_genes/{sample}.log",
    threads: 24
    shell:
        """
        featureCounts -s {params.strand} -a "{input.gtf}" {params.paired_end} \
            -T {threads} {params.extra} -o "{output}" "{input.bam}" \
            1> "{log}" 2>&1
        """
