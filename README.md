# BeringSea_shotgun_sequencing
Scripts for the compositional analysis of metagenomic shotgun sequencing data calssified by kraken2

Here, I provide the code that was used for the processing, filtering, and analysis of the metagenomic shotgun sequencing data of sedimentary ancient DNA from a marine sediment core from the Bering Sea for the paper by Zimmermann et al. 2023. 

The processing of the sequencing data and taxonomic classification were carried out on a HPC 
The code to process the data is part of this README. Further data processing, filtering and the analysis of the data was carried out in `R`. The R-scripts are part of this repository and the README guides through the order of running them.



Overview of used programs and versions
|Program|Version|Link|
|-------|------|------|
|FastQC|0.11.9|[Link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|FastUniq|1.1|[Link](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0052249)|
|Fastp|0.20.0|[Link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129281/)|
|kraken2|2.0.8-beta|[Link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)|
|HOPS|0.3.4|[Link](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1903-0#:~:text=HOPS%20is%20a%20versatile%20tool%20for%20high-throughput%20screening,enables%20large-scale%20metagenomic%20analyses%20of%20complex%20biological%20systems.)|


---


## Quality check of sequencing data

```
for i in *fastq.gz; do fastqc $i -t 4 -o fastqc_reports/; done
```

## Duplicate removal


## Data processing
Trimming low quality ends and removing reads with low complexity, removal of residual adapters and merging of overlapping forward and reverse reads.

```
fastp --in1=forward.fastq --in2=reverse.fastq --out1=unmerged_forward.fastq --out2=unmerged_reverse.fastq \
--unpaired1=unpaired_forward.fastq --unpaired2=unpaired_reverse.fastq --merge \
--merged_out=merged_samplename.fastq --failed_out=failed_samplename.fastq --thread=36 --verbose \
--adapter_fasta=adapter_2020_04_06.fasta --length_required=30 --overlap_len_require=10 --overlap_diff_percent_limit=20 \
--correction --overlap_diff_limit=5 --low_complexity_filter --complexity_threshold=30 --qualified_quality_phred=15 \
--unqualified_percent_limit=40 --n_base_limit=5 --json=samplename.json --html=samplename.html \
--report_title="samplename_fastp_report" --overrepresentation_analysis  --cut_front --cut_tail --cut_window_size=4 \
--cut_mean_quality=10 -q 10 -x 10
```

## Taxonomic classification with kraken2

[kraken2 manual](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)

### single-end mode for merged reads:
```
kraken2-2.0.8-beta/kraken2 --conf=0.2 --db kraken2_db merged_samplename.fastq --threads 36 \
--output samplename_merged.kraken --report samplename_merged.report
```

### paired-end mode for unmerged read pairs: 
```
kraken2-2.0.8-beta/kraken2 --conf=0.2 --paired --db kraken2_db unmerged_forward.fastq \
unmerged_reverse.fastq --threads 36 --report samplename_unmerged.report --output samplename_unmerged.kraken 
```


```
kraken2-inspect --threads 36 --db /home/ollie/projects/bio/db/kraken2/nt_2021_04_db  \
--use-mpa-style >> ~/krakeninspectdb/kraken2-inspect_nt_2021_04_db.txt
```

### Damage pattern analysis
