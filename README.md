

## Can-Seq script for identifying nucleotide transitions

### Installation

Requires:
- Linux or Mac OSX (likely Bash on Ubuntu on Windows too, but not tested)
- Anaconda / [Miniconda](https://conda.io/miniconda.html) installed
- [Bioconda channels added](http://ddocent.com//bioconda/)

Install the software in the can_seq environment:

```
conda create -n can_seq python=3.6 biopython trim-galore bowtie2 samtools varscan snakemake
```

- Download the can_seq script (green button above and right)

### Usage

1. Copy reference files (genomic region of each candidate gene) in genbank format to the ```gb``` directory.  Files must have the extension ```.gb```
2. Copy paired-end FASTQ sequences files to the ```raw``` directory.  These should be uncompressed and have the file extensions ```_1.fq``` and ```_2.fq``` for each of the pairs.  Example file/directory structure:
```
project-root----can_seq.snakefile
	     |--gb----gene_1.gb
	     |     |--gene_2.gb
	     |
	     |--raw----reads_1.fq
	            |--reads_2.fq 
```    
3. Activate the can_seq environment
```
conda activate can_seq
```
4. Run the script.  It should be copied to the root directory of the project
```
snakemake -s can_seq.snakefile --configfile config.json
```
### Essential requirements:

Genbank files covering the complete genomic region of the target gene should be correctly formatted
- To construct the protein sequence and determine amino acid substitutions, at least one set of CDS coordinates is required in the FEATURES section.
- e.g. ```     CDS             join(420..3134,3557..4432)```
- Genbank files sourced via NCBI-Entrez Gene gene should be correctly formatted (e.g. [RDR6](https://www.ncbi.nlm.nih.gov/nuccore/NC_003074.8?report=genbank&from=18348974&to=18353673&strand=true))

### What will happen

The script will trim adapters from the read files, align them to each reference, determine the location of any substitutions, generate CSV files with the results, and annotate new genbank (.gb) files with the location, substitution and frequency of the nucleotide substitution, and the position and amino acid substitution if any takes place.

Example CSV output for RDR6:
| CDS | Var nucl pos | Ref nucl |  Var nucl | Var percent | FWD Ref coverage | FWD Ref coverage | RVS Ref coverage | FWD Var coverage | FWD Rvs coverage | Var prot pos | Ref aa | Var aa |
|--|--|--|--|--|--|--|--|--|--|--|--|
1 | 475	| G	| A	| 5.49% | 3404 | 1982 | 192 | 121 | 19 | G | E |
1 | 1100 | G | A | 4.43% | 3276 | 2327 | 148 | 112 | 227 | W | * |
1 | 1545 | C | T | 3.59% | 3440 | 2261 | 131 | 81 | 376	R | * |
1 | 2474 | G | A | 3.41% | 3621 | 2661 | 146 | 76 | 685	W | * |
1 | 3427 | G | A | 5.29% | 3781 | 1999 | 205 | 118 | Intron | - | - |
