# BeringSea_shotgun_sequencing
Scripts for the compositional analysis of metagenomic shotgun sequencing data calssified by kraken2

Here, I provide the code that was used for the processing, filtering, and analysis of the metagenomic shotgun sequencing data of sedimentary ancient DNA from a marine sediment core from the Bering Sea for the paper by Zimmermann et al. 2023. 

The processing of the sequencing data and taxonomic classification were carried out on a HPC 
The code to process the data is part of this README. Further data processing, filtering and the analysis of the data was carried out in `R`. The R-scripts are part of this repository and the README guides through the order of running them.



Overview of used programs, versions used for the analysis, and the link to the publication
|Program|Version|Link to publication|
|-------|------|------|
|FastQC|0.11.9|[Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|FastUniq|1.1|[Link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249)|
|Fastp|0.20.0|[Link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/)|
|kraken2|2.0.8-beta|[Link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)|
|HOPS|0.3.4|[Link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1903-0#:~:text=HOPS%20is%20a%20versatile%20tool%20for%20high-throughput%20screening,enables%20large-scale%20metagenomic%20analyses%20of%20complex%20biological%20systems.)|
|R|4.0.3||


---


## Quality check of sequencing data

```
for i in *fastq.gz; do fastqc $i -t 4 -o fastqc_reports/; done
```

## Duplicate removal
```
fastuniq -i samplename.list -t q -o samplename_dedup_R1.fastq -p samplename_dedup_R2.fastq
```

## Data processing
Trimming low quality ends and removing reads with low complexity, removal of residual adapters and merging of overlapping forward and reverse reads.
Here is the link to the [fastp manual]() for further information.
```
fastp --in1=samplename_dedup_R1.fastq --in2=samplename_dedup_R2.fastq --out1=unmerged_forward.fastq --out2=unmerged_reverse.fastq \
--unpaired1=unpaired_forward.fastq --unpaired2=unpaired_reverse.fastq --merge \
--merged_out=merged_samplename.fastq --failed_out=failed_samplename.fastq --thread=36 --verbose \
--adapter_fasta=adapter_2020_04_06.fasta --length_required=30 --overlap_len_require=10 --overlap_diff_percent_limit=20 \
--correction --overlap_diff_limit=5 --low_complexity_filter --complexity_threshold=30 --qualified_quality_phred=15 \
--unqualified_percent_limit=40 --n_base_limit=5 --json=samplename.json --html=samplename.html \
--report_title="samplename_fastp_report" --overrepresentation_analysis  --cut_front --cut_tail --cut_window_size=4 \
--cut_mean_quality=10 -q 10 -x 10
```

## Taxonomic classification with kraken2

Here is the link to the [kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown) for further information.

### single-end mode for merged reads:
```
kraken2-2.0.8-beta/kraken2 --conf=0.2 --db kraken2_db merged_samplename.fastq --threads 36 \
--output samplename_merged.kraken --report samplename_merged.report
```

### paired-end mode for unmerged read pairs: 
```
kraken2-2.0.8-beta/kraken2 --conf=0.2 --paired --db kraken2_db unmerged_forward.fastq \
unmerged_reverse.fastq --threads 36 --report samplename_unmerged.report \
--output samplename_unmerged.kraken 
```

### Inspect the database for signal validation of the dominant fish families in our analysis
To check whether correlations could be influenced by a reference data bias of closely related taxa, we inspected the our custom full-nt kraken database.
Formal analysis was carried out afterwards using the R-script `fish_db_inspection.R`.
```
kraken2-inspect --threads 36 --db /home/ollie/projects/bio/db/kraken2/nt_2021_04_db  \
--use-mpa-style >> ~/krakeninspectdb/kraken2-inspect_nt_2021_04_db.txt
```

### Damage pattern analysis
Here is the link to the [HOPS manual]() for further information.
We prepared 3 datasets in which samples were binned according to their age to estimate whether damage patterns increase with time: set1 (), set2 (), set3 ().

```
#!/bin/bash 

#SBATCH --job-name=hops_nt_APMG689
#SBATCH -p xfat
#SBATCH --qos=large
#SBATCH -t 96:00:00 
#SBATCH --cpus-per-task=28
#SBATCH --mail-type=ALL
#SBATCH --mem=1400G 
#SBATCH --error="hopsNT-%j.err"
#SBATCH --out="hopsNT-%j.out"

#####Author: Lars Harms

# set variables (requires modification)
#================================================
WORKDIR=/shotgun/APMG689/hops/
INPUT=/shotgun/APMG689/output/out.fastp/191011_SND405_A_L007_APMG-6-1_fastp_merged_R2.fq.gz
OUTPUT=/shotgun/APMG689/hops/output/APMG-6-1
CONFIG=/shotgun/APMG689/hops/configfile.txt
INDEXDB=/db/hops/nt-hops-20-11-03_step8/
TAXFILE=/shotgun/APMG689/hops/taxalist_2022_05_06.txt
NCBIRESC=/db/hops/ncbi
MEM=1400

# preparing the working environment
#================================================
cd $WORKDIR
module load bio/hops/0.34

# create Hops config file
#================================================
echo "preProcess=0 
alignment=1 
dupRemOff=0 
filter=def_anc 
index=${INDEXDB}
pathToList=${TAXFILE}
resources=${NCBIRESC}

useSlurm = 0
threadsMalt=${SLURM_CPUS_PER_TASK}
maxMemoryMalt=${MEM}

threadsMaltEx=${SLURM_CPUS_PER_TASK}
maxMemoryMaltEx=${MEM}

threadsPost=${SLURM_CPUS_PER_TASK}
maxMemoryPost=256" > ${CONFIG}

# tasks to be performed
#================================================

srun hops -Xmx${MEM}G -input ${INPUT} -output ${OUTPUT} -m full -c ${CONFIG} 

```

---
# Data processing, filtering, compositional analysis and network analysis

The following R-scripts and python-script are run in the order they are listed below. Input files for running the scripts are included in the folder Input_files.

1. Cleaning up the kraken2-output and filtering for family level with  `1_prep_famlevel.R`
2. Making taxonomic subgroups based on manually checked taxa lists `2_prep_taxalists.R`
3. Resampling of reads to account for differences in per sample read counts for pelagic taxa `3_resampling_plankton_pelagic.R` and benthic taxa `3_resampling_benthic.R`. This script is based on the github script by Stefan Kruse [R-rarefaction](https://github.com/StefanKruse/R_Rarefaction).
4. Compositional analysis, Spearman correlations, and network analysis for pelagic taxa `4_analysis_pelagic.R` and benthic taxa `4_analysis_benthic.R`.
5. ``
6. ``
7. ``
8. ``



