# bifrost_assemblatron

The input data to this initial component will already be adapter - and quality trimmed data after "min_read_check", so this component when given a sample id and requirements checks will perform de novo assembly and a quality control step calculating various assembly metrics. 

## de novo Assembly (see pipeline.smk)
quality trimmed reads are stored within sample['categories']['paired_reads']['summary']['trimmed'] - so indexing is required to access then as done in the pipeline. With a basic spades command as such
```
spades -1 trimmed_read1.fastq.gz -2 trimmed_read1.fastq.gz --threads xx --isolate -o outputdirectory
```

As described in the documentation (https://ablab.github.io/spades/running.html) the "--isolate" is an assembly running mode: "This flag is highly recommended for high-coverage isolate and multi-cell Illumina data; improves the assembly quality and running time"

## rule assembly_qc
Takes the created de-novo assembly and filters the fasta files and calculates basic metrics for all of the split data files for future QC steps. The filtering steps follows:
1. Removes all contigs below length of 500 stored in (option --min-contig-len & --failed-length)
2. Removes all the contigs with a coverage below 10x (extracted from the spades header, option --cov-threshold & --failed-coverage)
3. Stores the passed files as the final filtered assembly files used for future components (--passed)

The command takes prefixes for the filtered files to generate both the filtered fasta files and statistical measurements
```
python rule__assembly_qc.py --assembly generated_de_novo_assembly.fasta --cov-threshold 10 --min-contig-len 500 --passed passed_file_prefix --failed-length filtered_on_length_prefix --failed-coverage filtered_on_coverage_prefix
```

# Rerun this sole component as a module for dev/test
```
snakemake --nolock --cores all xxx
```


