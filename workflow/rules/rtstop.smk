if RTSTOP_LIBRARY_SCHEMA is not None:

    rule rtstop_count:
        input:
            bam="results/align/{sample}.transcriptome.strand_filtered.unsorted.bam",
        output:
            "results/rtstop/{sample}.tsv.gz",
        log:
            "logs/rtstop/{sample}.log",
        conda:
            "../envs/rtstop.yaml"
        params:
            schema=RTSTOP_LIBRARY_SCHEMA,
        threads: 4
        script:
            "../scripts/rtstop_count.py"
