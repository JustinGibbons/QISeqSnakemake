#!/usr/bin/env bash

## getSequencerFromFastq.sh
# requires: fcid, gzcat

## run this script in the directory with your fastq.gz files (or specify input directory with -d option) to determine the sequencing machine based on the header line from each fastq file. The script prints one line per FASTQ file, with tab-separated fields: sampleName    flowcell_id    sequencing_machine


set -euo pipefail

usage() {
  cat <<'EOF'
Usage: getSequencerFromFastq.sh [-d DIR]

Options:
  -d, --dir   Directory containing *.fastq.gz files (default: current directory)
  -h, --help  Show this help and exit

Output:
  Prints tab-separated: sampleName  flowcell_id  sequencing_machine
EOF
}

# If no arguments were provided, display help and exit
if [[ $# -eq 0 ]]; then
  usage
  exit 0
fi

indir="."
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d|--dir)
      [[ $# -ge 2 ]] || { echo "Error: Missing argument for $1" >&2; exit 2; }
      indir="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Error: Unknown argument: $1" >&2; usage; exit 2;;
  esac
done

# Hardcode gzcat
ZCAT="gzcat"

# Ensure fcid is available
command -v fcid >/dev/null 2>&1 || { echo "Error: 'fcid' not found in PATH." >&2; exit 1; }

# Validate directory
[[ -d "$indir" ]] || { echo "Error: Directory not found: $indir" >&2; exit 1; }

shopt -s nullglob

found_any=false
for fq in "$indir"/*.fastq.gz; do
  [[ -e "$fq" ]] || continue
  found_any=true

  filename="$(basename -- "$fq")"
  sampleName="${filename%.fastq.gz}"

  # Read exactly the first line without causing SIGPIPE
  header=""
  if ! IFS= read -r header < <("$ZCAT" "$fq" 2>/dev/null); then
    echo "Warning: Could not read header from '$fq'." >&2
    continue
  fi

  if [[ -z "$header" || "${header:0:1}" != "@" ]]; then
    echo "Warning: '$fq' does not look like a valid FASTQ (missing '@' header)." >&2
    continue
  fi

  # Extract flowcell ID (3rd colon-delimited field)
  # Example: '@ERR14906234.1 LH00275:27:222VT7LT1:2:2207:35257:8386/2' -> '222VT7LT1'
  flowcell_id="$(printf '%s' "$header" | cut -d':' -f3)"

  if [[ -z "$flowcell_id" ]]; then
    echo "Warning: Could not parse flowcell ID from '$fq'." >&2
    continue
  fi

  # Run fcid (stdin or arg)
  if machine="$(printf '%s\n' "$flowcell_id" | fcid 2>/dev/null)"; then
    :
  elif machine="$(fcid "$flowcell_id" 2>/dev/null)"; then
    :
  else
    echo "Warning: 'fcid' failed for '$flowcell_id' from '$fq'." >&2
    continue
  fi

  printf "%s\t%s\t%s\n" "$sampleName" "$flowcell_id" "$machine"
done

if ! $found_any; then
  echo "No *.fastq.gz files found in: $indir" >&2
  exit 1
fi
