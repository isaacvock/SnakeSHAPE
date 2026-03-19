# Then can make coverage tracks
rule coverage_tracks:
    input:
        bam="results/filter/{sample}.bam",
        bai="results/filter/{sample}.bam.bai",
    output:
        "results/coverage/{sample}.bw",
    log:
        "logs/coverage_tracks/{sample}.log",
    conda:
        "../envs/genomictools.yaml"
    params:
        extra=BAMCOVERAGE_EXTRA,
    threads: 8
    shell:
        """
        bamCoverage -p {threads} -b "{input.bam}" -o "{output}" {params.extra} \
            1> "{log}" 2>&1
        """
