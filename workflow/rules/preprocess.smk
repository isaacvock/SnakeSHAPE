if TRIMMER == "cutadapt":

    rule cutadapt_pe:
        input:
            get_units_fastqs,
        output:
            fastq1="results/trimmed/{sample}_R1.fastq.gz",
            fastq2="results/trimmed/{sample}_R2.fastq.gz",
            qc="results/trimmed/{sample}.qc.txt",
        params:
            adapters=config.get("cutadapt_adapters", ""),
            extra=config.get("cutadapt_extra", ""),
        log:
            "logs/cutadapt_pe/{sample}.log",
        threads: 8
        wrapper:
            "v7.5.0/bio/cutadapt/pe"

    rule cutadapt_se:
        input:
            get_units_fastqs,
        output:
            fastq="results/trimmed/{sample}_single.fastq.gz",
            qc="results/trimmed/{sample}.qc.txt",
        params:
            adapters=config.get("cutadapt_adapters", ""),
            extra=config.get("cutadapt_extra", ""),
        log:
            "logs/cutadapt_se/{sample}.log",
        threads: 8
        wrapper:
            "v7.5.0/bio/cutadapt/se"

elif TRIMMER == "fastp":

    # Single-end data
    rule fastp_se:
        input:
            sample=get_units_fastqs,
        output:
            trimmed="results/trimmed/{sample}_single.fastq.gz",
            failed="results/trimmed/{sample}_single.failed.fastq.gz",
            html="results/trimmed/{sample}_single.html",
            json="results/trimmed/{sample}.json",
        log:
            "logs/fastp_se/{sample}.log",
        params:
            adapters=config.get("fastp_adapters"),
            extra=config.get("fastp_extra"),
        threads: 8
        wrapper:
            "v7.2.0/bio/fastp"

    # Paired-end data
    rule fastp_pe:
        input:
            sample=get_units_fastqs,
        output:
            trimmed=[
                "results/trimmed/{sample}_R1.fastq.gz",
                "results/trimmed/{sample}_R2.fastq.gz",
            ],
            # Unpaired reads separately
            unpaired1="results/trimmed/{sample}.unpaired.R1.fastq.gz",
            unpaired2="results/trimmed/{sample}.unpaired.R2.fastq.gz",
            failed="results/trimmed/{sample}_paired.failed.fastq.gz",
            html="results/trimmed/{sample}.html",
            json="results/trimmed/{sample}.json",
        log:
            "logs/fastp_pe/{sample}.log",
        params:
            adapters=config.get("fastp_adapters"),
            extra=config.get("fastp_extra"),
        threads: 8
        wrapper:
            "v7.2.0/bio/fastp"

else:
    raise ValueError("trimmer should be cutadapt or fastp!")


rule fastqc:
    input:
        unpack(get_fq),
    output:
        qc=directory("results/fastqc/{sample}"),
    log:
        "logs/fastqc/{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    params:
        extra=FASTQC_EXTRA,
    resources:
        mem_mb=9000,
    threads: 4
    shell:
        """
        mkdir -p {output.qc}
        fastqc {params.extra} --threads {threads} --outdir "{output.qc}" {input} \
            1> "{log}" 2>&1
        """
