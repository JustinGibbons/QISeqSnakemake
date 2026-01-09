#!/usr/bin/env python
"""
tradis_fixbases_normalized_correctDup.py

Python translation of the Perl script 'tradis.fixbases_normalized_correctDup.pl'.
This utility consolidates insertion counts by proximity, computes a "mapped"
total based on a set of pre-approved ("OK") insertion sites, and emits normalized
counts per 100,000 mapped reads for anchor positions that are also OK.

Behavior mirrors the original Perl logic:

Inputs
------
1) Required CLI argument: a path to the "OK sites" file. Each line must be:
   "<amount> <chrom> <pos>"
   For every 'pos' in this file, the window (pos-30 .. pos+30) is marked as OK.
2) STDIN: lines of "<amount> <chrom> <pos>" from your aggregation pipeline.

Aggregation & merging
---------------------
- Uses a 1,000 bp bin gate per chromosome: the first position observed in a bin
  becomes its anchor and is kept in output order.
- Subsequent positions landing in the same bin are merged into an existing
  nearby position within ±1..±10 bp (preferring +1..+10, then -1..-10).
- Counts are accumulated at the resolved anchor/nearby position.

Mapped total & output
---------------------
- 'mapped' is the sum of amounts that fall within the OK-site windows.
- Prints 'mapped' on its own line.
- Then, for each anchor in first-seen order, if the anchor position itself is OK:
  prints: "<normalized>\t<count>\t<chrom>\t<pos>"
  where normalized = round((count / mapped) * 100000, 2).

Notes
-----
- If 'mapped' is 0 (no OK hits), normalization is reported as 0.00 to avoid
  division-by-zero (this is a safety enhancement; Perl would error out).

Author: (Python translation by Copilot)
"""

import sys
import re
from collections import defaultdict
from typing import Dict, Set, Tuple, List


def load_ok_sites(path: str) -> Dict[str, Set[int]]:
    """
    Load OK sites from file and expand each pos to a +/-30 bp window.

    Each input line: "<amount> <chrom> <pos>" (amount is ignored; presence suffices).
    For every line, positions {pos-30 .. pos+30} are marked as OK.

    Parameters
    ----------
    path : str
        Path to the OK sites file.

    Returns
    -------
    Dict[str, Set[int]]
        Mapping from chromosome -> set of OK positions (expanded windows).
    """
    ok: Dict[str, Set[int]] = defaultdict(set)
    try:
        with open(path, "r") as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                # Strip leading whitespace like Perl's s/^\s+//g
                line = re.sub(r"^\s+", "", line)
                if not line:
                    continue
                parts = re.split(r"\s+", line)
                if len(parts) < 3:
                    continue
                # amount is not used; only chrom and pos matter for OK marking
                chrom = parts[1]
                try:
                    pos = int(parts[2])
                except ValueError:
                    continue
                for delta in range(-30, 31):
                    ok[chrom].add(pos + delta)
    except OSError as e:
        sys.stderr.write(f"ERROR: cannot open OK sites file '{path}': {e}\n")
        sys.exit(1)

    return ok


def find_pos(pos: int, chrom: str, counts: Dict[str, Dict[int, int]]) -> int:
    """
    Find an existing nearby position within ±1..±10 bp to merge into.

    Preference order:
      +1, +2, ... +10, then -1, -2, ... -10
    Returns 'pos' itself if no nearby existing position is found.

    Parameters
    ----------
    pos : int
        Candidate position to merge.
    chrom : str
        Chromosome identifier.
    counts : Dict[str, Dict[int, int]]
        Accumulated counts: counts[chrom][pos] -> total amount.

    Returns
    -------
    int
        Nearby existing position to merge into, or the original 'pos'.
    """
    for i in range(1, 11):
        if (pos + i) in counts[chrom]:
            return pos + i
        if (pos - i) in counts[chrom]:
            return pos - i
    return pos


def main():
    # --- Parse required OK sites argument (mirrors Perl: die if missing) ---
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python tradis_fixbases_normalized_correctDup.py <ok_sites_file>\n")
        sys.exit(1)

    ok_sites_path = sys.argv[1]
    ok_sites = load_ok_sites(ok_sites_path)

    # Data structures mirroring Perl:
    counts: Dict[str, Dict[int, int]] = defaultdict(lambda: defaultdict(int))
    # around[chrom][pos // 1000] marks that a 1,000 bp bin has an anchor already
    around: Dict[str, Dict[int, bool]] = defaultdict(dict)
    # Maintain first-seen anchor order as a list of (chrom, pos)
    order: List[Tuple[str, int]] = []

    mapped_total = 0

    # --- Read all stdin lines (Perl pushes into @readin, then loops) ---
    stdin_lines = [ln.rstrip("\n") for ln in sys.stdin]

    # --- Process lines: strip leading whitespace, parse amount chrom pos ---
    for line in stdin_lines:
        line = re.sub(r"^\s+", "", line)
        if not line:
            continue
        parts = re.split(r"\s+", line)
        if len(parts) < 3:
            continue

        try:
            amount = int(parts[0])
            chrom = parts[1]
            pos = int(parts[2])
        except ValueError:
            # Skip malformed lines
            continue

        pos_bin = pos // 1000

        if pos_bin in around[chrom]:
            # Bin already has an anchor: merge into nearby existing position
            pos_new = find_pos(pos, chrom, counts)
            counts[chrom][pos_new] += amount
            # Add to mapped_total only if the resolved pos_new is OK
            if pos_new in ok_sites.get(chrom, set()):
                mapped_total += amount
        else:
            # First time we touch this bin: create anchor and remember order
            counts[chrom][pos] = counts[chrom][pos] + amount
            around[chrom][pos_bin] = True
            order.append((chrom, pos))
            # Add to mapped_total if this anchor is OK
            if pos in ok_sites.get(chrom, set()):
                mapped_total += amount

    # --- Output mapped_total (Perl prints it on its own line) ---
    print(f"{mapped_total}")

    # --- Output anchors in first-seen order, but only if anchor is OK ---
    for chrom, pos in order:
        if pos in ok_sites.get(chrom, set()):
            total = counts[chrom][pos]
            if mapped_total > 0:
                norm = round((total / mapped_total) * 100000, 2)
            else:
                # Safety enhancement vs. Perl: avoid division-by-zero
                norm = 0.00
            # Match Perl: normalized (2 decimals), then raw count, chrom, pos
            print(f"{norm:.2f}\t{total}\t{chrom}\t{pos}")


if __name__ == "__main__":
    main()

