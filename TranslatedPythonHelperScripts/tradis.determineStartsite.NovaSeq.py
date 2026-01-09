#!/usr/bin/env python
"""
Determine transposon insertion start sites from SAM input.

Behavior mirrors the Perl script tradis.determineStartsite.NovaSeq.pl:

- Reads SAM-format lines from stdin.
- Uses CIGAR string + strand flag to calculate insertion position.
- Uses an offset that must be provided as a command‑line argument.
- For reverse‐strand mappings:
    Length = M + D − I (excluding leading soft‑clipped bases)
    Final position = pos + length + offset
- For forward‑strand mappings:
    Final position = pos − offset

Outputs:
    <chrom> <calculated_start> <length_or_0> <cigar>

Author: (Python translation by Copilot)
"""

import sys
import re


def get_length(cigar: str) -> int:
    """
    Calculate effective aligned length from a CIGAR string
    following the Perl script's rules:

    Rules replicated:
    - Strip initial soft clip (e.g., "12S")
    - M (match) adds length
    - D (deletion) adds length
    - I (insertion) SUBTRACTS length
    - Subsequent soft-clips are ignored

    Parameters
    ----------
    cigar : str
        A standard SAM CIGAR string.

    Returns
    -------
    int
        Effective mapped length.
    """
    # Remove initial soft‐clip (e.g. "14S")
    cigar = re.sub(r"^\d+S", "", cigar)

    length = 0

    # Parse step-by-step until CIGAR fully consumed
    while cigar:
        m = re.match(r"^(\d+)M", cigar)
        if m:
            length += int(m.group(1))
            cigar = cigar[m.end():]
            continue

        m = re.match(r"^(\d+)I", cigar)
        if m:
            length -= int(m.group(1))
            cigar = cigar[m.end():]
            continue

        m = re.match(r"^(\d+)D", cigar)
        if m:
            length += int(m.group(1))
            cigar = cigar[m.end():]
            continue

        m = re.match(r"^(\d+)S", cigar)
        if m:
            cigar = cigar[m.end():]
            continue

        # If stuck, avoid infinite loop
        break

    return length


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python determine_startsite.py <offset>\n")
        sys.exit(1)

    offset = int(sys.argv[1])

    for line in sys.stdin:
        line = line.rstrip()
        if not line or line.startswith("@"):
            continue  # skip headers or empty

        fields = line.split("\t")
        if len(fields) < 6:
            print(line[:80])
            continue

        flag = int(fields[1])
        chrom = fields[2]
        pos = int(fields[3])
        cigar = fields[5]

        # Reverse strand if flag bit 0x10 is set
        reverse_strand = bool(flag & 0x010)

        # Case 1: reverse strand, CIGAR ends with M
        if reverse_strand and re.search(r"\d+M$", cigar):
            lgth = get_length(cigar)
            start = pos + lgth + offset
            print(f"{chrom}\t{start}\t{lgth}\t{cigar}")

        # Case 2: forward strand match
        elif re.search(r"\d+M", cigar):
            start = pos - offset
            print(f"{chrom}\t{start}\t0")

        else:
            # Fallback — match Perl behavior
            print(line[:80])


if __name__ == "__main__":
    main()

