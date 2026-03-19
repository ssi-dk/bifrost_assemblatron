# bifrost_assemblatron

This component is used to perform a *de novo* assembly from trimmed paired-end sequence reads obtained from another [component](https://github.com/ssi-dk/bifrost_min_read_check/tree/master). Following the *de novo* assembly, basic filtering is performed and metrics is calculated, to ensure high quality contigs.  

## Requirements
- The component creates the *de novo* assembly using the tool [spades](https://github.com/ablab/spades).
- The versions are described in the [environment.yaml](https://github.com/ssi-dk/bifrost_assemblatron/blob/master/environment.yml)
- The alignment uses the trimmed reads from another [component](https://github.com/ssi-dk/bifrost_min_read_check/tree/master).

## Download
```bash
git clone https://github.com/ssi-dk/bifrost_assemblatron.git
cd bifrost_assemblatron
git submodule init
git submodule update
bash install.sh -i LOCAL
conda activate bifrost_assemblatron_vx.x.x
export BIFROST_INSTALL_DIR='/your/path/'
BIFROST_DB_KEY="/your/key/here/" python -m bifrost_assemblatron -h
```
## Run the snakemake analysis
Each component can be run on each sample individually using one snakemake command, replacing the string passed to the **--config sample_name=" "** with the appropriate dataset name. The provided **component_name=** takes as an argument *<component_name>__<version_number>*. The component name aligns with the GitHub repo name, which is structured like *bifrost_<component_name>* (e.g. *bifrost_assemblatron* -> component name *assemblatron*), and the version number aligns with the current [GitHub tag](https://github.com/ssi-dk/bifrost_assemblatron/tags) / or conda environment [version](https://github.com/ssi-dk/bifrost_assemblatron/blob/master/setup.py) (e.g. *v.0.0.2*) defined during the bifrost component setup. 
```bash
snakemake -p --nolock --cores 5 -s <github_path>/pipeline.smk --config sample_name="insert sample name" component_name=assemblatron__v2.3.3 --rerun-incomplete
```

## Analysis
### do novo assembly
The spades command is defined in the [pipeline](https://github.com/ssi-dk/bifrost_assemblatron/blob/master/bifrost_assemblatron/pipeline.smk), using the assembly running mode *"--isolate"* ([spades documentation](https://ablab.github.io/spades/running.html)) is recommended for high-coverage isolate and multi-cell Illumina data, and this improves the assembly quality and decreases the running time.

One example of the assembled contig headers
```bash
>24000006_MW-ESCEC_1_length_505018_cov_24.431370
AACCTGCGACCAATTGATTAAAAGTCAACTGCTCTACCAACTGAGCTAACGACCCACTTT
TTCGTTGCTTTCGGTTTGTTTGATATCCCGTGGCAACGGCGGCATATATTACTGATTTCA
GACTTGAGCGCAACAAAAATTTCGATGTAGATCACTCAACTGCTTATGATTCGCACGACA
```
From the contig header, the depth of coverage can be extracted and used for QC on the assembled contigs. 

### rule assembly_qc
The following filtering and metric calculations are performed in the [rule__assembly_qc](https://github.com/ssi-dk/bifrost_assemblatron/blob/master/bifrost_assemblatron/rule__assembly_qc.py) using the following steps:
1. Removes all contigs below the length of 500 bp (option --min-contig-len & --failed-length)
2. Removes all the contigs with a coverage below 10x (extracted from the spades header, option --cov-threshold & --failed-coverage)
3. Stores the passed files as the final filtered assembly files used for future components (--passed)

The command takes prefixes for the filtered files to generate both the filtered fasta files and statistical measurements
```
python rule__assembly_qc.py --assembly generated_de_novo_assembly.fasta --cov-threshold 10 --min-contig-len 500 --passed passed_file_prefix --failed-length filtered_on_length_prefix --failed-coverage filtered_on_coverage_prefix
```

The calculated metrics include:
- shortest contig length
- total number of filtered contigs
- total length of contigs
- N50 : which is a measurement representing the shortest contig length required to cover 50% of the total assembly length
- Average depth of coverage from the filtered contigs

With one example of an output statistical file shown below. 
```
group   cov_threshold   min_contig_len  no_contigs      sum_len n50     mean_cov
cov_ge_10       10      500     95      5018390 186344  24.88
```


