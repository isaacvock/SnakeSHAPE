if config.get("aligner", "star") == "star":

    # Build STAR index
    rule index:
        input:
            fasta=config["genome"],
            gtf=config["annotation"],
        output:
            outdir=directory(config.get("indices", "star_index")),
            outfile=star_indices,
        params:
            extra=config["star_index_params"],
        log:
            "logs/index/star_index_genome.log",
        conda:
            "../envs/star.yaml"
        threads: 24
        script:
            "../scripts/star-index.py"

    # Align with STAR
    # I think passing GTF here is redundant given how
    # indexing is done, but I don't think it hurts and
    # can catch provided-index-edge-cases
    rule align:
        input:
            unpack(get_fq),
            idx=config.get("indices", "star_index"),
            gtf=config["annotation"],
        output:
            aln="results/align/{sample}.bam",
            transcriptome_raw="results/align/{sample}.transcriptome.unsorted.bam",
            reads_per_gene="results/align/{sample}_ReadsPerGene.out.tab",
            log="results/align/{sample}_Log.out",
            sj="results/align/{sample}_SJ.out.tab",
            log_final="results/align/{sample}_Log.final.out",
        log:
            "logs/align/{sample}_star.log",
        conda:
            "../envs/star.yaml"
        params:
            out_prefix=lambda wc: f"results/align/{wc.sample}_",
            reads=get_cli_input_fastqs,
            read_cmd=get_gzip_read_command,
            extra=config.get("star_align_params", ""),
        threads: 24
        shell:
            """
            STAR --runThreadN {threads} \
                --genomeDir "{input.idx}" \
                --readFilesIn {params.reads} \
                {params.read_cmd} \
                --outFileNamePrefix "{params.out_prefix}" \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode GeneCounts TranscriptomeSAM \
                --quantTranscriptomeBan Singleend \
                --sjdbGTFfile "{input.gtf}" \
                {params.extra} \
                1> "{log}" 2>&1

            mv "{params.out_prefix}Aligned.sortedByCoord.out.bam" "{output.aln}"
            mv "{params.out_prefix}Aligned.toTranscriptome.out.bam" "{output.transcriptome_raw}"
            mv "{params.out_prefix}ReadsPerGene.out.tab" "{output.reads_per_gene}"
            mv "{params.out_prefix}Log.out" "{output.log}"
            mv "{params.out_prefix}SJ.out.tab" "{output.sj}"
            mv "{params.out_prefix}Log.final.out" "{output.log_final}"
            """

    rule sort_transcriptome_bam:
        input:
            bam="results/align/{sample}.transcriptome.unsorted.bam",
        output:
            "results/align/{sample}.transcriptome.bam",
        log:
            "logs/align/{sample}_transcriptome_sort.log",
        conda:
            "../envs/genomictools.yaml"
        threads: 8
        shell:
            """
            samtools sort -@ {threads} -o "{output}" "{input.bam}" \
                1> "{log}" 2>&1
            """

    rule rsem_prepare_reference:
        input:
            fasta=config["genome"],
            gtf=config["annotation"],
        output:
            grp=RSEM_REFERENCE_GRP,
            transcriptome=RSEM_TRANSCRIPTOME_FASTA,
        log:
            "logs/rsem/reference.log",
        conda:
            "../envs/rsem.yaml"
        params:
            prefix=RSEM_REFERENCE_PREFIX,
            extra=RSEM_PREPARE_REFERENCE_EXTRA,
        threads: 24
        shell:
            """
            mkdir -p "$(dirname "{params.prefix}")"
            rsem-prepare-reference --gtf "{input.gtf}" --star -p {threads} \
                {params.extra} "{input.fasta}" "{params.prefix}" \
                1> "{log}" 2>&1
            """

    rule rsem_quantify:
        input:
            unpack(get_fq),
            ref=RSEM_REFERENCE_GRP,
        output:
            genes="results/rsem/{sample}.genes.results",
            isoforms="results/rsem/{sample}.isoforms.results",
        log:
            "logs/rsem/{sample}.log",
        conda:
            "../envs/rsem.yaml"
        params:
            prefix=RSEM_REFERENCE_PREFIX,
            paired_end=lambda wc: "--paired-end" if is_paired_end(wc.sample) else "",
            strandedness=RSEM_STRANDEDNESS,
            reads=get_cli_input_fastqs,
            gzip_flag=get_rsem_gzip_flag,
            extra=RSEM_QUANTIFY_EXTRA,
            sample_prefix=lambda wc: f"results/rsem/{wc.sample}",
        threads: 24
        shell:
            """
            mkdir -p "$(dirname "{params.sample_prefix}")"
            rsem-calculate-expression --star --no-bam-output {params.paired_end} \
                --strandedness {params.strandedness} {params.gzip_flag} \
                -p {threads} {params.extra} {params.reads} "{params.prefix}" \
                "{params.sample_prefix}" \
                1> "{log}" 2>&1
            """

### Filtering out unmapped, secondary, QC-failed, duplicate, and supplementary reads.
### For paired-end data, also drop reads with unmapped mates.
rule filter:
    input:
        "results/align/{sample}.bam",
    output:
        "results/filter/{sample}.bam",
    log:
        "logs/filter/{sample}.log",
    conda:
        "../envs/genomictools.yaml"
    params:
        flags=lambda wc: get_filter_samtools_flags(wc.sample),
    threads: 8
    shell:
        """
        samtools view -@ {threads} -bh {params.flags} -o "{output}" "{input}" \
            1> "{log}" 2>&1
        """
