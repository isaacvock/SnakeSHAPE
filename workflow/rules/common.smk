import glob
import shlex

import pandas as pd
from snakemake.utils import validate

import os
from pathlib import Path
from snakemake.exceptions import WorkflowError

##### CHECK VALIDITIY OF SAMPLES CSV #####


### Path to samples
SAMPLES_PATH = Path(config.get("samples"))


### Does samples file exist?
if not SAMPLES_PATH.exists():
    raise WorkflowError(f"Samples CSV not found!")


### Load sample information
samples = (
    pd.read_csv(SAMPLES_PATH, dtype={"sample_name": str})
    .set_index(["sample_name"], drop=False)
    .sort_index()
)


##### HELPFUL CONFIG PARSING #####

# Which adapter trimmer to use
TRIMMER = config.get("trimmer", "cutadapt")

# featureCounts settings
FEATURECOUNTS_EXTRA = config.get(
    "featurecounts_extra",
    "-O --nonOverlap 0 --extraAttributes gene_id",
)
FASTQC_EXTRA = config.get("fastqc_extra", "")
PERBASE_EXTRA = config.get("perbase_extra", "")
RSEM_REFERENCE_PREFIX = config.get(
    "rsem_reference_prefix",
    "resources/rsem/reference/rsem",
)
RSEM_REFERENCE_GRP = f"{RSEM_REFERENCE_PREFIX}.grp"
RSEM_TRANSCRIPTOME_FASTA = f"{RSEM_REFERENCE_PREFIX}.transcripts.fa"
RSEM_TRANSCRIPTOME_FASTA_FAI = f"{RSEM_TRANSCRIPTOME_FASTA}.fai"
RSEM_PREPARE_REFERENCE_EXTRA = config.get("rsem_prepare_reference_extra", "")
RSEM_QUANTIFY_EXTRA = config.get("rsem_quantify_extra", "")

# Need some sort of tangible file to confirm that alignment
# index has been built/exists
if config.get("indices", "index").endswith("/"):
    INDEX_PATH = str(config.get("indices", "index"))
    INDEX_PATH = INDEX_PATH[:-1]
else:
    INDEX_PATH = str(config.get("indices", "index"))
star_indices = glob.glob(f"{INDEX_PATH}/genomeParameters.txt")


def has_config_value(value):
    return pd.notna(value) and str(value).strip() != ""


def validate_bamcoverage_extra(extra):
    tokens = shlex.split(str(extra))
    forbidden = {"--extendReads", "-e"}

    if any(token in forbidden for token in tokens):
        raise WorkflowError(
            "bamCoverage read extension is not appropriate for spliced RNA-seq "
            "coverage tracks. Remove --extendReads or -e from bamcoverage_extra."
        )

    return extra


def parse_bool(value):
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False

    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n", ""}:
        return False

    raise WorkflowError(f"Could not interpret boolean config/sample value: {value}")


def get_featurecounts_strand():
    strandedness = str(config.get("strandedness", "reverse")).strip().lower()
    strand_map = {
        "unstranded": 0,
        "none": 0,
        "forward": 1,
        "sense": 1,
        "reverse": 2,
        "antisense": 2,
    }

    if strandedness not in strand_map:
        valid = ", ".join(sorted(strand_map))
        raise WorkflowError(
            f"Unsupported strandedness '{strandedness}'. Expected one of: {valid}"
        )

    return strand_map[strandedness]


def get_rsem_strandedness():
    strandedness = str(config.get("strandedness", "reverse")).strip().lower()
    strand_map = {
        "unstranded": "none",
        "none": "none",
        "forward": "forward",
        "sense": "forward",
        "reverse": "reverse",
        "antisense": "reverse",
    }

    if strandedness not in strand_map:
        valid = ", ".join(sorted(strand_map))
        raise WorkflowError(
            f"Unsupported strandedness '{strandedness}'. Expected one of: {valid}"
        )

    return strand_map[strandedness]


def get_rtstop_library_schema():
    schema = str(config.get("rtstop_library_schema", "")).strip().lower()
    if schema == "":
        return None

    aliases = {
        "single_end": "single_end",
        "single": "single_end",
        "se": "single_end",
        "read1": "read1",
        "r1": "read1",
        "read2": "read2",
        "r2": "read2",
    }

    if schema not in aliases:
        valid = ", ".join(sorted(aliases))
        raise WorkflowError(
            f"Unsupported rtstop_library_schema '{schema}'. Expected one of: {valid}"
        )

    return aliases[schema]


FEATURECOUNTS_STRAND = get_featurecounts_strand()
RSEM_STRANDEDNESS = get_rsem_strandedness()
RTSTOP_LIBRARY_SCHEMA = get_rtstop_library_schema()
BAMCOVERAGE_EXTRA = validate_bamcoverage_extra(
    config.get("bamcoverage_extra", config.get("deeptools_params", ""))
)

##### FUNCTIONS USED THROUGHOUT #####


### Get paths to fastqs for trimming
def get_units_fastqs(wildcards):
    s = samples.loc[(wildcards.sample)]
    if not has_config_value(s["fq1"]):
        # SRA sample
        accession = s["sra"]
        if is_paired_end(wildcards.sample):
            return expand(
                "sra/{accession}_{read}.fastq",
                accession=accession,
                read=["R1", "R2"],
            )
        else:
            return [f"sra/{accession}_R1.fastq"]
    if not is_paired_end(wildcards.sample):
        return [
            s["fq1"],
        ]
    else:
        return [s["fq1"], s["fq2"]]


### Figure out if sample is paired end
def is_paired_end(sample):
    sample_units = samples.loc[sample]

    if has_config_value(sample_units["fq2"]):
        return True

    if has_config_value(sample_units["sra"]):
        return parse_bool(sample_units.get("PE", False))

    return False


### Get paths to fastqs for alignment
def get_fq(wildcards):
    if config["trim_reads"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "results/trimmed/{sample}_{read}.fastq.gz",
                        read=["R1", "R2"],
                        **wildcards,
                    ),
                )
            )
        # single end sample
        return {"fq1": "results/trimmed/{sample}_single.fastq.gz".format(**wildcards)}
    else:
        # no trimming, use raw reads
        fqs = get_units_fastqs(wildcards)
        if len(fqs) == 1:
            return {"fq1": f"{fqs[0]}"}
        elif len(fqs) == 2:
            return {"fq1": f"{fqs[0]}", "fq2": f"{fqs[1]}"}
        else:
            raise ValueError(f"Expected one or two fastq file paths, but got: {fqs}")


def get_input_fastqs(wildcards):
    return list(get_fq(wildcards).values())


def get_cli_input_fastqs(wildcards):
    return " ".join(f'"{fq}"' for fq in get_input_fastqs(wildcards))


def get_gzip_read_command(wildcards):
    fastqs = get_input_fastqs(wildcards)
    gzipped = [str(fq).endswith(".gz") for fq in fastqs]

    if any(gzipped) and not all(gzipped):
        raise WorkflowError(
            f"Sample '{wildcards.sample}' mixes gzipped and plain-text FASTQs, which "
            "is not supported."
        )

    return "--readFilesCommand gunzip -c" if all(gzipped) else ""


def get_rsem_gzip_flag(wildcards):
    fastqs = get_input_fastqs(wildcards)
    gzipped = [str(fq).endswith(".gz") for fq in fastqs]

    if any(gzipped) and not all(gzipped):
        raise WorkflowError(
            f"Sample '{wildcards.sample}' mixes gzipped and plain-text FASTQs, which "
            "is not supported."
        )

    return "--star-gzipped-read-file" if all(gzipped) else ""


def get_filter_samtools_flags(sample):
    flags = [
        "-F 0x4",  # unmapped read
        "-F 0x100",  # secondary alignment
        "-F 0x200",  # vendor/platform QC failure
        "-F 0x400",  # duplicate
        "-F 0x800",  # supplementary alignment
    ]

    if is_paired_end(sample):
        flags.append("-F 0x8")  # mate unmapped

    return " ".join(flags)


def strip_suffix(path, suffix):
    path = str(path)
    if not path.endswith(suffix):
        raise WorkflowError(f"Expected path ending in '{suffix}', got '{path}'")
    return path[: -len(suffix)]


def get_rsem_reference_prefix_from_grp(path):
    return strip_suffix(path, ".grp")


def get_rsem_sample_prefix_from_genes(path):
    return strip_suffix(path, ".genes.results")


### What samples get merged
SAMPS_TO_MERGE = samples["sample_name"].tolist()


def validate_rtstop_library_schema():
    if RTSTOP_LIBRARY_SCHEMA is None:
        return

    paired_samples = [sample for sample in SAMPS_TO_MERGE if is_paired_end(sample)]
    single_samples = [sample for sample in SAMPS_TO_MERGE if not is_paired_end(sample)]

    if paired_samples and single_samples:
        raise WorkflowError(
            "RT-stop counting currently expects all samples to share one library "
            "layout. Mixed single-end and paired-end samples are not supported when "
            "rtstop_library_schema is set."
        )

    if paired_samples and RTSTOP_LIBRARY_SCHEMA not in {"read1", "read2"}:
        raise WorkflowError(
            "For paired-end RT-stop libraries, rtstop_library_schema must be "
            "'read1' or 'read2'."
        )

    if single_samples and RTSTOP_LIBRARY_SCHEMA != "single_end":
        raise WorkflowError(
            "For single-end RT-stop libraries, rtstop_library_schema must be "
            "'single_end'."
        )


validate_rtstop_library_schema()


### Get final output
def get_final_output():
    ### What needs to be included?
    # 1) FastQC reports
    # 2) Genome-aligned BAMs
    # 3) Transcriptome-aligned BAMs
    # 4) Coverage tracks
    # 5) featureCounts gene quantification
    # 6) RSEM transcript quantification
    # 7) perbase output
    # 8) RT-stop counts (when enabled)

    final_output = []
    final_output.extend(expand("results/fastqc/{sample}", sample=SAMPS_TO_MERGE))
    final_output.extend(expand("results/align/{sample}.bam", sample=SAMPS_TO_MERGE))
    final_output.extend(
        expand("results/align/{sample}.transcriptome.bam", sample=SAMPS_TO_MERGE)
    )
    final_output.extend(expand("results/coverage/{sample}.bw", sample=SAMPS_TO_MERGE))
    final_output.extend(
        expand(
            "results/quantify_genes/{sample}_gene.featureCounts",
            sample=SAMPS_TO_MERGE,
        )
    )
    final_output.extend(
        expand("results/rsem/{sample}.genes.results", sample=SAMPS_TO_MERGE)
    )
    final_output.extend(
        expand("results/rsem/{sample}.isoforms.results", sample=SAMPS_TO_MERGE)
    )
    final_output.extend(
        expand(
            "results/perbase/{sample}.tsv.gz",
            sample=SAMPS_TO_MERGE,
        )
    )
    if RTSTOP_LIBRARY_SCHEMA is not None:
        final_output.extend(
            expand(
                "results/rtstop/{sample}.tsv.gz",
                sample=SAMPS_TO_MERGE,
            )
        )

    return final_output
