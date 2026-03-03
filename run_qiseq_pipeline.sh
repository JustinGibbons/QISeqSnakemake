#!/bin/bash


# Define the usage function
usage() {
    echo "" >&2
    echo "Usage: bash $(basename "$0") [-n <integer>]" >&2
    echo "  -n <integer>: Number of jobs that can run in parallel. A common choice is the number of samples (required)" >&2
    echo "  -s Which snakefile to run. Typically either snakefile_forward or snakfile_reverse depending on library prep method, but can be your own version (required)" >&2
    echo "  -h          : Display this help message" >&2
    echo "" >&2
    echo "" >&2
    echo "    First edit this script to change the --mail-user value to YOUR email address (~line 66)" >&2
    echo "    Then run this script in the same directory as your updated snakefile." >&2
    echo "    The pipeline will run quietly in the background. Check smart_slurm.out and the ./Logs directory for progress." >&2
    echo "" >&2
    exit 1
}
    

## Parse command-line options
# Variables to hold command line arguments
SAMPLES=""
SNAKE_FILE=""
EMAIL=""
# Use a leading ':' to make getopts quiet, letting you handle errors
while getopts ":hn:s:" opt; do
  case "$opt" in
    h)
      usage
      ;;
    n)
      SAMPLES="$OPTARG"
      ;;
    s) 
      SNAKE_FILE="$OPTARG"
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


if [[ -z "$SNAKE_FILE" ]]; then
      echo "Error: snakefile  must be specified with the -s option." >&2
        usage
    fi


##### Run Snakemake with SLURM cluster submission #####
mkdir -p Logs && nohup snakemake -j $SAMPLES --keep-going --snakefile $SNAKE_FILE --cluster '
sbatch -p rra --qos rra \
-n 1 \
--cpus-per-task {threads} \
--mem {resources.mem_mb} \
--time {resources.runtime} \
--output Logs/{rule}.{wildcards}.o \
--error Logs/{rule}.{wildcards}.e \
--mail-type END,FAIL \
--mail-user jgibbons1@usf.edu \
--job-name {rule}.{wildcards}' --default-resources runtime=30 >qiseq_snakemake_stout.out 2>&1 &



