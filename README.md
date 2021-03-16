# bifrost_assemblatron

This component is run given a sample id already added into the bifrostDB. From this it'll pull the paired_reads and preform assembly of the data and remapping of the raw reads to the assembly. From the resulting files various QC metrics will be captured.

## Programs: (see Dockerfile) 
```
snakemake-minimal==5.31.1;
bbmap==38.58;
skesa==2.4.0;
minimap2==2.17;
samtools==1.11;
cyvcf2==0.30.1;
quast==5.0.2;
```

## Summary of c run: (see pipeline.smk and config.yaml)
```
java -ea -cp /opt/conda/opt/bbmap-38.58-0/current/ jgi.BBDuk in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo qtrim=r minlength=30 1> {log.out_file} 2> {log.err_file}
skesa --use_paired_ends --fastq {input.filtered_reads} --contigs_out {output.contigs} 1> {log.out_file} 2> {log.err_file}
quast.py {input.contigs} -o {output.quast} 1> {log.out_file} 2> {log.err_file}
bbsketch.sh in={input.contigs} out={output.sketch} 1> {log.out_file} 2> {log.err_file}
stats.sh {input.contigs} 1> {log.out_file} 2> {log.err_file}
minimap2 --MD -ax sr {input.contigs} {input.filtered_reads} 1> {output.mapped} 2> {log.err_file}
samtools stats {input.mapped} 1> {output.stats} 2> {log.err_file}
pileup.sh in={input.mapped} basecov={output.coverage} out={output.pileup} 1> {log.out_file} 2> {log.err_file}
rule__summarize_depth.py
callvariants.sh in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}
rule__summarize_variants.py
sed -e 's/Contig/{params.sample_name}/' {input.contigs} > {output.contigs}
```
