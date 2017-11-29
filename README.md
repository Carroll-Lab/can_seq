## Can-Seq script for identifying nucleotide transitions

### Installation

Requires:
- Linux or Mac OSX
- Anaconda / Miniconda installed
- Bioconda channels added

Install the software in the can_seq environment:

```
conda create -n can_seq python=3.6 biopython trim-galore bowtie2 samtools varscan snakemake
```

- Download the can_seq script

### Usage

1. Copy reference files (genomic region of each candidate gene) in genbank format to the ```gb``` directory.  Files must have the extension ```.gb```
2. Copy paired-end FASTQ sequences files to the ```raw``` directory.  These should be uncomplressed and have the file extensions ```_1.fq``` and ```_2.fq``` for each of the pairs
3. Run the script.  It should be in the root sirectory of the project
```
snakemake -s can_seq.snakefile
```

### What will happen

The script will trim adapters from the read files, align them to each reference, determine the location of any substitutions, generate CSV files with the results, and annotate new genbank (.gb) files with the location, substitution and frequency of the substitution


