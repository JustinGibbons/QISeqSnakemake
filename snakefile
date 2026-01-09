## BEFORE RUNNING: make sure you've activated the conda environment containing all QIseq pipeline dependencies first, like so:
##	module load hub.apps/anaconda/2022.12
##	conda activate QISeq_Deps_Snakemake7
import os
import glob

####################### USER INPUTS ###########################
# provide directory path to your fastq.gz files:
experiment_name="48917_L5678" 
SampleDir = "/shares/pi_ja2/JG/48917_L5678/FASTQs"
	## note that your fastq filenames must end in "_R1.fastq.gz" and "_R2.fastq.gz"

# provide directory path to your reference genome:
RefDir = "RefGenome"
	# RefDir = "/shares/omicshub/References/Plasmodium/Pfalciparum/nf54"

# provide complete path to reference genome fasta file:
refFasta = RefDir+"/PlasmoDB-56_PknowlesiA1H1_Genome.fasta"
	# refFasta = RefDir+"/PfNF54.fasta"

# set offset parameter depending on sequencing platform. Default is -2 for NovaSeq
# NovaSeq (both NovaSeqX and NovaSeq6000):
Offset=-2

## HiSeq:
#Offset=4

## MiSeq (old platform; not sure about newest one):
#Offset=5
###############################################################


#ScriptsDir = "/shares/pi_ja2/JG/NewQISeqWorkflows/Updated_QISeq_Helper_Scripts"
insertSize=1000 ##Thomas said he thought was 1,000
Min_Unique=1 ##Minimum number of unique insertions
merge_distance=3 ##Merge predicted insertion sites if within merge_distance bp


determineStartsite_exec="TranslatedPythonHelperScripts/tradis.determineStartsite.NovaSeq.py"
fixbases_exec="TranslatedPythonHelperScripts/tradis_fixbases.py"
fixbases_normalized_correctDup_exec="TranslatedPythonHelperScripts/tradis_fixbases_normalized_correctDup.py"
combine_tables_exec="TranslatedPythonHelperScripts/mergeToTable_copilot_refactored.py"
def get_ref_genome_name(inputpath):
    file_base_name = os.path.basename(inputpath)
    genome_name = ".".join(file_base_name.split(".")[0:-1])
    return genome_name


def get_sample_names(indir):
    fastqfiles=glob.glob(indir.rstrip("/")+"/*.fastq.gz")
    sample_names=[]

    for fastq in fastqfiles:
        fastq_basename=os.path.basename(fastq)
        sample_name=fastq_basename.replace(".fastq.gz","")
        sample_name="_".join(sample_name.split("_")[0:-1])
        sample_names.append(sample_name)

    return(sample_names)

genomeName = get_ref_genome_name(refFasta)

#samples = ["MP_IC0_F1_D10", "MP_IC10_F2_D10"]
samples=set(get_sample_names(SampleDir))

rule all:
    input:
        expand("FinalCounts/{sample}",sample=samples),
        expand("Insertion_Counts/{sample}.txt",sample=samples),
        expand("SamtoolStats/{sample}_stats.txt",sample=samples)

rule delete_adapter_from_readsF:
    priority: -1
    input:
        SampleDir+"/{sample}_R1.fastq.gz",
    output:
        temp("TrimmedR1/{sample}_R1_adapter_trimmed.fastq.gz"),
    threads: 20
    resources:
        runtime=60
    params:
        "TrimmedR1/{sample}_R1_adapter_trimmed.fastq",
    shell:
        """zcat {input} | perl -e 'while (<>){{ print; $_=<STDIN>; print substr($_,6);$_=<STDIN>;print; $_=<STDIN>; print substr($_,6); }}' > {params} && pigz -p {threads} {params}"""

rule create_bowtie_index:
    priority: 1
    input:
        refFasta

    output:
        RefDir+"/"+genomeName+".1.bt2"

    shell:
        "bowtie2-build {input} {RefDir}/{genomeName}"

rule samtools_index_genome:
    priority: 1
    input:
        refFasta

    output:
        refFasta+".fai"

    shell:
        "samtools faidx {input}"

rule bowtie_align:
    priority: -1
    input:
        R1="TrimmedR1/{sample}_R1_adapter_trimmed.fastq.gz",
        R2=SampleDir+"/{sample}_R2.fastq.gz",
        ref_genome= RefDir+"/"+genomeName+".1.bt2"

    output:
        temp("BowtieSams/{sample}.sam")

    threads: 20

    resources:
        runtime=120

    shell:
        "bowtie2 -X {insertSize} -p {threads} --very-sensitive -N 1 -L 31 --rdg 5,2 -x {RefDir}/{genomeName} -1 {input.R1} -2 {input.R2} -S {output}"


rule samTobam:
    priority: 1
    input:
        samfile="BowtieSams/{sample}.sam",
        genomeIndex=refFasta+".fai"
    output:
        temp("BowtieBams/{sample}.bam")

    shell:
        """awk -va={wildcards.sample} '{{if ($1~/^@/) {{print}} else  {{print $0"\tRG:Z:"a}}}}' {input.samfile} | samtools view -b -t {input.genomeIndex} - > {output}"""


rule create_new_bam_headers:
    priority: 1
    ##Not sure why this is neccessary.
    ##Also something to do with read groups. May be able to combine
    ##the previous rule and the subsquent 1 using an script or Picard command
    input:
        "BowtieBams/{sample}.bam"
    output:
        temp("BamHeaders/{sample}_bam_headers.txt")
    shell:
        """samtools view -H {input} >{output}
        echo -e "@RG\tID:{wildcards.sample}\tSM:1" >>{output}"""

rule add_new_bam_headers:
    priority: 1
    input:
        bam="BowtieBams/{sample}.bam",
        header="BamHeaders/{sample}_bam_headers.txt"

    output:
        ##Can mark as temporary
        temp("New_Header_Bams/{sample}.bam")

    shell:
        "samtools reheader {input.header} {input.bam} >{output}"


rule sort_bam_file:
    priority: 1
    input:
        "New_Header_Bams/{sample}.bam"
        ##Can mark as temporarty
    output:
        temp("Sorted_Bams/{sample}.bam")

    resources:
        runtime=60

    shell:
        "samtools sort {input} -o {output}"


rule mark_duplicates:
    priority: 5
    input:
        "Sorted_Bams/{sample}.bam"
        ##These are not temporary
    output:
        "MarkDups_Bams/{sample}.bam"
    resources:
        runtime=60,
        mem_mb=lambda wildcards, input: max(2*input.size_mb,4000)
	## JO upped max mem from 2000 to 4000

    shell:
        "picard MarkDuplicates VALIDATION_STRINGENCY=SILENT M=/dev/null ASSUME_SORTED=TRUE MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I={input} O={output}"


rule index_bams:
    priority: 1
    input:
        "MarkDups_Bams/{sample}.bam"

    output:
        "MarkDups_Bams/{sample}.bai"

    shell:
        "samtools index {input} {output}"


rule count_insertions_exclude_duplicates:
    priority: 1
    input:
        bam="MarkDups_Bams/{sample}.bam",
        ##This is just to make sure index gets made
        index="MarkDups_Bams/{sample}.bai"
    output:
        "Insertion_Counts_No_Dups/{sample}.txt"

    shell:
        """samtools view -f 66 -F 1024 {input.bam} | awk ' !($6 ~ "S")' | {determineStartsite_exec} {Offset} | sort | uniq -c |sort -rn | {fixbases_exec} {Min_Unique} > {output}"""


rule count_insertions:
    priority: 10
    input:
        bam="MarkDups_Bams/{sample}.bam",
        ##This is just to make sure index gets made
        index="MarkDups_Bams/{sample}.bai",
        unique_count="Insertion_Counts_No_Dups/{sample}.txt"

    output:
        "Insertion_Counts/{sample}.txt"

    shell:
        """samtools view -f 66 {input.bam} | {determineStartsite_exec} {Offset} |  sort | uniq -c |sort -rn | {fixbases_normalized_correctDup_exec}  {input.unique_count}  > {output}"""

rule run_samtools_stats:
    priority: 1
    input:
        bam="MarkDups_Bams/{sample}.bam",
        index="MarkDups_Bams/{sample}.bai"

    output:
        "SamtoolStats/{sample}_stats.txt"
    threads: 10
    resources:
        runtime=10,
        mem_mb=1000
    shell:
        "samtools stats -@ {threads} {input.bam} > {output}"


rule combine_count_tables:
    priority: 1
    
    input:
        expand("Insertion_Counts/{sample}.txt",sample=samples)

    output:
        temp(touch(expand("FinalCounts/{sample}",sample=samples)))

    params:
        "Insertion_Counts"

    shell:
        "{combine_tables_exec} {params} {merge_distance} > FinalCounts/{experiment_name}_insertion_counts.tsv"
