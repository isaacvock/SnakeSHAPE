rule get_sra_se:
    wildcard_constraints:
        accession=r"[^_/]+",
    output:
        "sra/{accession}.fastq",
    log:
        "logs/get_sra/{accession}.log",
    conda:
        "../envs/sra.yaml"
    threads: 6
    shell:
        """
        mkdir -p sra "$(dirname "{log}")"
        fasterq-dump --skip-technical -e {threads} --outdir sra "{wildcards.accession}" \
            1> "{log}" 2>&1
        """


rule get_sra_pe:
    wildcard_constraints:
        accession=r"[^_/]+",
    output:
        fastq1="sra/{accession}_1.fastq",
        fastq2="sra/{accession}_2.fastq",
    log:
        "logs/get_sra/{accession}.log",
    conda:
        "../envs/sra.yaml"
    threads: 6
    shell:
        """
        mkdir -p sra "$(dirname "{log}")"
        fasterq-dump --skip-technical -e {threads} --outdir sra "{wildcards.accession}" \
            1> "{log}" 2>&1

        orphan="sra/{wildcards.accession}.fastq"
        if [ -e "$orphan" ]; then
            rm -f "$orphan"
        fi
        """
