from collections import Counter
import gzip
from pathlib import Path

import pysam


def get_rtstop_site(read):
    if read.is_reverse:
        return read.reference_end
    return read.reference_start + 1


def informative_read(read, schema):
    if schema == "single_end":
        return True
    if schema == "read1":
        return read.is_read1
    if schema == "read2":
        return read.is_read2
    raise ValueError(f"Unsupported RT-stop schema: {schema}")


schema = snakemake.params.schema
bam_path = Path(snakemake.input.bam)
output_path = Path(snakemake.output[0])
log_path = Path(snakemake.log[0])

counts = Counter()
total_alignments = 0
counted_events = 0

with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam:
        total_alignments += 1

        if read.is_unmapped or read.is_supplementary or read.is_qcfail:
            continue

        if not informative_read(read, schema):
            continue

        site = get_rtstop_site(read)
        counts[(read.reference_name, site)] += 1
        counted_events += 1

output_path.parent.mkdir(parents=True, exist_ok=True)
with gzip.open(output_path, "wt") as handle:
    handle.write("transcript_id\tposition_1based\trtstop_count\n")
    for (transcript_id, position), count in sorted(counts.items()):
        handle.write(f"{transcript_id}\t{position}\t{count}\n")

log_path.parent.mkdir(parents=True, exist_ok=True)
with log_path.open("w") as handle:
    handle.write(f"schema\t{schema}\n")
    handle.write(f"input_bam\t{bam_path}\n")
    handle.write(f"total_alignments\t{total_alignments}\n")
    handle.write(f"counted_events\t{counted_events}\n")
    handle.write(f"unique_sites\t{len(counts)}\n")
