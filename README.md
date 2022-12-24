# BeringSea_shotgun_sequencing
Scripts for the compositional analysis of metagenomic shotgun sequencing data calssified by kraken2

Here, I provide the code that was used for the processing, filtering, and analysis of the metagenomic shotgun sequencing data of sedimentary ancient DNA from a marine sediment core from the Bering Sea for the paper by Zimmermann et al. 2023. 

The processing of the sequencing data and taxonomic classification were carried out on a HPC 


Overview of used programs and versions
|Program|Version|
|-------|------|
|FastQC|0.11.9|
|FastUniq|1.1|
|Fastp|0.20.0|
|kraken2|2.0.8-beta|
|HOPS|0.3.4|


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
