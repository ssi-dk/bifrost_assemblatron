# bifrost_assemblatron

This component is run given a sample id already added into the bifrostDB. From this it'll pull the paired_reads and contigs and remap the raw reads to the assembly. From the resulting files various QC metrics will be captured.

## Programs: (see Dockerfile) 
```
snakemake-minimal==5.31.1;
bbmap==38.58;
skesa==2.4.0;
```

## Summary of c run: (see pipeline.smk and config.yaml)
```
java -ea -cp /opt/conda/opt/bbmap-38.58-0/current/ jgi.BBDuk in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo qtrim=r minlength=30 1> {log.out_file} 2> {log.err_file}
skesa --use_paired_ends --fastq {input.filtered_reads} --contigs_out {output.contigs} 1> {log.out_file} 2> {log.err_file}
```
