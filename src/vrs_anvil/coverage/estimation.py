"""Synthetic region profiles inferred from variant positions, not read coverage."""

import gzip
from collections import defaultdict
from pathlib import Path


def vcf_to_bed(vcf_path: str | Path, bed_path: str | Path, *, window: int = 0) -> None:
    """Write the merged union of VCF POS windows as headerless BED3.

    POS=100 with window=10 becomes [89, 110). All FILTER states and
    variant ALT types are included, irrespective of genotypes or INFO/END.
    Placeholder-only ALT records are excluded. This profile does not establish
    sequencing coverage or callability. Input intervals are sorted in memory.

    Raises:
        ValueError: For invalid windows, VCF structure, or identical paths.
        OSError: For file access errors.

    """
    if isinstance(window, bool) or not isinstance(window, int) or window < 0:
        msg = "window must be a nonnegative integer"
        raise ValueError(msg)
    source, destination = Path(vcf_path), Path(bed_path)
    if source.resolve() == destination.resolve() or (
        destination.exists() and source.samefile(destination)
    ):
        msg = "Input and output paths must differ"
        raise ValueError(msg)
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with source.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    header: list[str] | None = None
    with opener(source, "rt") as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.rstrip("\r\n").split("\t")
            error = f"Invalid VCF at line {line_number}"
            if line_number == 1:
                if not line.startswith("##fileformat=VCFv"):
                    raise ValueError(error + ": missing VCF fileformat header")
                continue
            if line.startswith("##") and header is None:
                continue
            if line.startswith("#CHROM") and header is None:
                expected = [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                ]
                if fields[:8] != expected or (
                    len(fields) > len(expected)
                    and (len(fields) < len(expected) + 2 or fields[8] != "FORMAT")
                ):
                    raise ValueError(error + ": invalid column header")
                header = fields
                continue
            if header is None or len(fields) != len(header) or line.startswith("#"):
                raise ValueError(error + ": missing header or incorrect column count")
            if not fields[0] or any(not field for field in fields[:8]):
                raise ValueError(error + ": empty required field")
            if not fields[1].isascii() or not fields[1].isdigit() or int(fields[1]) < 1:
                raise ValueError(error + ": POS must be a positive integer")
            if all(
                alt in {".", "*", "<NON_REF>", "<*>"} for alt in fields[4].split(",")
            ):
                continue
            position = int(fields[1])
            intervals[fields[0]].append(
                (max(0, position - 1 - window), position + window)
            )
    if header is None:
        msg = "Invalid VCF: missing column header"
        raise ValueError(msg)
    with destination.open("w") as output:
        for chrom in sorted(intervals):
            start, end = -1, -1
            for next_start, next_end in sorted(intervals[chrom]):
                if start >= 0 and next_start > end:
                    output.write(f"{chrom}\t{start}\t{end}\n")
                    start, end = next_start, next_end
                else:
                    start = next_start if start < 0 else start
                    end = max(end, next_end)
            output.write(f"{chrom}\t{start}\t{end}\n")
