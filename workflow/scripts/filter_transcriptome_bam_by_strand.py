from pathlib import Path

import pysam


def parse_gtf_attributes(field):
    attributes = {}
    for raw_item in field.strip().split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if " " not in item:
            continue
        key, value = item.split(" ", 1)
        attributes[key] = value.strip().strip('"')
    return attributes


def normalize_strandedness(strandedness):
    normalized = str(strandedness).strip().lower()
    mapping = {
        "forward": "forward",
        "sense": "forward",
        "reverse": "reverse",
        "antisense": "reverse",
        "none": "unstranded",
        "unstranded": "unstranded",
    }

    if normalized not in mapping:
        valid = ", ".join(sorted(mapping))
        raise ValueError(
            f"Unsupported strandedness '{strandedness}'. Expected one of: {valid}"
        )

    return mapping[normalized]


def get_read_key(read):
    if read.is_paired:
        if read.is_read1:
            return (read.query_name, "R1")
        if read.is_read2:
            return (read.query_name, "R2")
        raise ValueError(
            f"Paired read '{read.query_name}' is missing read1/read2 flag information."
        )
    return (read.query_name, "SE")


def opposite_strand(strand):
    return "-" if strand == "+" else "+"


def infer_transcript_strand(read, strandedness):
    if strandedness == "unstranded":
        return None

    aligned_strand = "-" if read.is_reverse else "+"

    if read.is_paired:
        if read.is_read1:
            same_as_transcript = strandedness == "forward"
        elif read.is_read2:
            same_as_transcript = strandedness == "reverse"
        else:
            raise ValueError(
                f"Paired read '{read.query_name}' is missing read1/read2 flags."
            )
    else:
        same_as_transcript = strandedness == "forward"

    return aligned_strand if same_as_transcript else opposite_strand(aligned_strand)


def load_transcript_strands(annotation_path):
    transcript_strands = {}

    with open(annotation_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            strand = fields[6]
            if strand not in {"+", "-"}:
                continue

            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if transcript_id is None:
                continue

            previous = transcript_strands.get(transcript_id)
            if previous is not None and previous != strand:
                raise ValueError(
                    f"Transcript '{transcript_id}' has conflicting strands "
                    f"('{previous}' vs '{strand}') in annotation."
                )

            transcript_strands[transcript_id] = strand

    return transcript_strands


strandedness = normalize_strandedness(snakemake.params.strandedness)
annotation_path = Path(snakemake.input.annotation)
genome_bam_path = Path(snakemake.input.genome_bam)
transcriptome_bam_path = Path(snakemake.input.transcriptome_bam)
output_bam_path = Path(snakemake.output[0])
log_path = Path(snakemake.log[0])

transcript_strands = load_transcript_strands(annotation_path)
genome_read_strands = {}

genome_records = 0
with pysam.AlignmentFile(genome_bam_path, "rb") as genome_bam:
    for read in genome_bam:
        genome_records += 1
        genome_read_strands[get_read_key(read)] = infer_transcript_strand(
            read, strandedness
        )

total_transcriptome_records = 0
kept_records = 0
dropped_missing_genome = 0
dropped_strand_mismatch = 0

output_bam_path.parent.mkdir(parents=True, exist_ok=True)
with pysam.AlignmentFile(transcriptome_bam_path, "rb") as transcriptome_bam, pysam.AlignmentFile(
    output_bam_path, "wb", template=transcriptome_bam
) as output_bam:
    for read in transcriptome_bam:
        total_transcriptome_records += 1

        read_key = get_read_key(read)
        inferred_strand = genome_read_strands.get(read_key)
        if read_key not in genome_read_strands:
            dropped_missing_genome += 1
            continue

        if inferred_strand is None:
            output_bam.write(read)
            kept_records += 1
            continue

        transcript_id = read.reference_name
        transcript_strand = transcript_strands.get(transcript_id)
        if transcript_strand is None:
            raise ValueError(
                f"Transcript '{transcript_id}' from transcriptome BAM was not found "
                "in the annotation transcript_id field."
            )

        if transcript_strand != inferred_strand:
            dropped_strand_mismatch += 1
            continue

        output_bam.write(read)
        kept_records += 1

log_path.parent.mkdir(parents=True, exist_ok=True)
with log_path.open("w", encoding="utf-8") as handle:
    handle.write(f"annotation\t{annotation_path}\n")
    handle.write(f"genome_bam\t{genome_bam_path}\n")
    handle.write(f"transcriptome_bam\t{transcriptome_bam_path}\n")
    handle.write(f"strandedness\t{strandedness}\n")
    handle.write(f"genome_records\t{genome_records}\n")
    handle.write(f"genome_read_keys\t{len(genome_read_strands)}\n")
    handle.write(f"transcriptome_records\t{total_transcriptome_records}\n")
    handle.write(f"kept_records\t{kept_records}\n")
    handle.write(f"dropped_missing_genome\t{dropped_missing_genome}\n")
    handle.write(f"dropped_strand_mismatch\t{dropped_strand_mismatch}\n")
