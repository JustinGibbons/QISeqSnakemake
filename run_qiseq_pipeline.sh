#!/bin/bash


# Define the usage function
usage() {
    echo "" >&2
    echo "Usage: bash $(basename "$0") [-n <integer>]" >&2
    echo "  -n <integer>: Number of samples to process with the QIseq pipeline (required)" >&2
    echo "  -h          : Display this help message" >&2
    echo "" >&2
    echo "" >&2
    echo "    First edit this script to change the --mail-user value to YOUR email address (~line 66)" >&2
    echo "    Then run this script in the same directory as your updated snakefile." >&2
    echo "    The pipeline will run quietly in the background. Check smart_slurm.out and the ./Logs directory for progress." >&2
    echo "" >&2
    exit 1
}
# INSTRUCTIONS: change the --mail-user value on line 66 to YOUR email address!!
    
    ## note you must also have updated 'snakefile' in the current directory with your user inputs for the pipeline to run.

## Parse command-line options
# Variable to hold the integer argument n
SAMPLES=""

# Use a leading ':' to make getopts quiet, letting you handle errors
while getopts ":hn:" opt; do
  case "$opt" in
    h)
      usage
      ;;
    n)
      SAMPLES="$OPTARG"
      ;;
    :) 
      echo "Error: Option -$OPTARG requires an argument." >&2
      usage
      ;;
    \?) 
      echo "Error: Invalid option -$OPTARG" >&2
      usage
      ;;
  esac
done

# Shift to remove parsed options and their arguments from the positional parameters
shift $((OPTIND-1))

# Check if the required argument was provided
if [[ -z "$SAMPLES" ]]; then
  echo "Error: Number of samples must be specified with the -n option." >&2
  usage
fi


##### Run Snakemake with SLURM cluster submission #####
mkdir -p Logs && nohup snakemake -j $SAMPLES --keep-going --snakefile snakefile --cluster '
sbatch -p rra --qos rra \
-n 1 \
--cpus-per-task {threads} \
--mem {resources.mem_mb} \
--time {resources.runtime} \
--output Logs/{rule}.{wildcards}.o \
--error Logs/{rule}.{wildcards}.e \
--mail-type END,FAIL \
--mail-user jgibbons1@usf.edu \
--job-name {rule}.{wildcards}' --default-resources runtime=30 >smart_slurm2.out 2>&1 &



