#!/usr/bin/env python
"""
tradis_fixbases.py

Python translation of the Perl script 'tradis.fixbases.pl' used in QISeq-style
workflows to consolidate insertion site counts by proximity. The logic mirrors
the original behavior:

- Input: lines of "<amount> <chrom> <pos>" from stdin (usually from aggregation
  like: cut -f1-3 | sort | uniq -c | sort -rn).
- Binning: For each chromosome, tracks whether a 1,000 bp bin (pos // 1000)
  has been encountered. If the bin is new, a position is stored directly.
  If the bin was seen earlier, attempts to merge the current count into an
  existing nearby position within ±1..±10 bp on the same chromosome.
- Merging rule: Prefer the first match found at +1..+10 bp; if none, try
  -1..-10 bp. If a nearby position exists, add the count to that position;
  otherwise, keep the original position.
- Output: Prints "<amount>\t<chrom>\t<pos>" for all positions whose total
  counts are ≥ minAmount (default 1), sorted by chromosome then position.

Author: (Python translation by Copilot)
"""

import sys
import re
from collections import defaultdict
from typing import Dict, Tuple


def find_pos(pos: int, chrom: str, counts: Dict[str, Dict[int, int]]) -> int:
    """
    Return an existing nearby position to which the current count should be merged.
    Search order: +1..+10, then -1..-10. If none found, return the original pos.

    Parameters
    ----------
    pos : int
        The candidate position to merge.
    chrom : str
        Chromosome identifier.
    counts : dict[str, dict[int, int]]
        Nested dict of accumulated counts: counts[chrom][pos] = amount.

    Returns
    -------
    int
        Either a nearby existing position within ±10 bp or the original pos.
    """
    # Try positive offsets first, then negative offsets
    for i in range(1, 11):
        if pos + i in counts[chrom]:
            return pos + i
        if pos - i in counts[chrom]:
            return pos - i
    return pos


def main():
    # Parse optional minAmount from argv; default is 1
    # (Matches Perl default behavior when arg is missing.)
    if len(sys.argv) >= 2:
        try:
            min_amount = int(sys.argv[1])
        except ValueError:
            sys.stderr.write("ERROR: minAmount must be an integer.\n")
            sys.exit(1)
    else:
        min_amount = 1

    # Data structures (mirror Perl hashes):
    # counts[chrom][pos] -> total amount
    counts: Dict[str, Dict[int, int]] = defaultdict(lambda: defaultdict(int))

    # around[chrom][pos // 1000] -> seen bin marker (value = True)
    around: Dict[str, Dict[int, bool]] = defaultdict(dict)

    # Read from stdin
    for raw in sys.stdin:
        line = raw.rstrip("\n")

        # Strip leading whitespace (Perl's s/^\s+//g)
        line = re.sub(r"^\s+", "", line)
        if not line:
            continue

        # Expect three fields: amount, chrom, pos (split on whitespace)
        parts = re.split(r"\s+", line)
        if len(parts) < 3:
            # Be forgiving—skip malformed lines
            continue

        try:
            amount = int(parts[0])
            chrom = parts[1]
            pos = int(parts[2])
        except ValueError:
            # Skip lines that don't parse cleanly
            continue

        pos_bin = pos // 1000

        if pos_bin in around[chrom]:
            # We've seen this 1,000 bp bin for this chromosome:
            # attempt to merge into an existing nearby position (±10 bp).
            pos_new = find_pos(pos, chrom, counts)
            counts[chrom][pos_new] += amount
        else:
            # First occurrence of this bin: store this exact position
            counts[chrom][pos] = counts[chrom][pos] + amount
            around[chrom][pos_bin] = True

    # Emit positions with total counts >= min_amount, sorted by chrom then pos
    for chrom in sorted(counts.keys()):
        for pos in sorted(counts[chrom].keys()):
            total = counts[chrom][pos]
            if total >= min_amount:
                print(f"{total}\t{chrom}\t{pos}")


if __name__ == "__main__":
    main()

