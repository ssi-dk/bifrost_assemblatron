#- Templated section: start ------------------------------------------------------------------------
import os
import sys
import traceback

from bifrostlib import common
from bifrostlib.datahandling import SampleReference
from bifrostlib.datahandling import Sample
from bifrostlib.datahandling import ComponentReference
from bifrostlib.datahandling import Component
from bifrostlib.datahandling import SampleComponentReference
from bifrostlib.datahandling import SampleComponent
os.umask(0o2)

try:
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample:Sample = Sample.load(sample_ref) # schema 2.1
    if sample is None:
        raise Exception("invalid sample passed")
    component_ref = ComponentReference(name=config['component_name'])
    component:Component = Component.load(reference=component_ref) # schema 2.1
    if component is None:
        raise Exception("invalid component passed")
    samplecomponent_ref = SampleComponentReference(name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference()))
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    if samplecomponent is None:
        samplecomponent:SampleComponent = SampleComponent(sample_reference=sample.to_reference(), component_reference=component.to_reference()) # schema 2.1
    common.set_status_and_save(sample, samplecomponent, "Running")
except Exception as error:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not samplecomponent.has_requirements():
        common.set_status_and_save(sample, samplecomponent, "Requirements not met")
    if samplecomponent['status'] == "Running":
        common.set_status_and_save(sample, samplecomponent, "Failure")

rule all:
    input:
        # file is defined by datadump function
        f"{component['name']}/datadump_complete"
    run:
        common.set_status_and_save(sample, samplecomponent, "Success")

rule setup:
    output:
        init_file = touch(temp(f"{component['name']}/initialized")),
    params:
        folder = component['name']
    run:
        samplecomponent['path'] = os.path.join(os.getcwd(), component['name'])
        samplecomponent.save()


rule_name = "check_requirements"
rule check_requirements:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        folder = rules.setup.output.init_file,
    output:
        check_file = f"{component['name']}/requirements_met",
    params:
        samplecomponent
    run:
        if samplecomponent.has_requirements():
            with open(output.check_file, "w") as fh:
                fh.write("")

#- Templated section: end --------------------------------------------------------------------------

#* Dynamic section: start **************************************************************************
rule_name = "setup__filter_reads_with_bbduk"
rule setup__filter_reads_with_bbduk:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.check_requirements.output.check_file,
        reads = sample['categories']['paired_reads']['summary']['data']
    output:
        filtered_reads = temp(f"{component['name']}/filtered.fastq")
    params:
        adapters = component['resources']['adapters_fasta']  # This is now done to the root of the continuum container
    shell:
        "java -ea -cp /opt/conda/opt/bbmap-38.58-0/current/ jgi.BBDuk in={input.reads[0]} in2={input.reads[1]} out={output.filtered_reads} ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tbo qtrim=r minlength=30 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly__skesa"
rule assembly__skesa:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        filtered_reads = rules.setup__filter_reads_with_bbduk.output.filtered_reads,
    output:
        contigs = f"{component['name']}/contigs.fasta"
    shell:
        "skesa --use_paired_ends --fastq {input.filtered_reads} --contigs_out {output.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__quast_on_contigs"
rule assembly_check__quast_on_contigs:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__skesa.output
    output:
        quast = directory(f"{component['name']}/quast")
    shell:
        "quast.py {input.contigs} -o {output.quast} 1> {log.out_file} 2> {log.err_file}"


rule_name = "assembly_check__sketch_on_contigs"
rule assembly_check__sketch_on_contigs:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__skesa.output
    output:
        sketch = f"{component['name']}/contigs.sketch"
    shell:
        "bbsketch.sh in={input.contigs} out={output.sketch} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__stats"
rule post_assembly__stats:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    message:
        "Running step: {rule}"
    input:
        contigs = rules.assembly__skesa.output
    output:
        stats = touch(f"{component['name']}/post_assermbly__stats")
    shell:
        "stats.sh {input.contigs} 1> {log.out_file} 2> {log.err_file}"


rule_name = "post_assembly__mapping"
rule post_assembly__mapping:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__skesa.output,
        filtered_reads = rules.setup__filter_reads_with_bbduk.output.filtered_reads
    output:
        mapped = temp(f"{component['name']}/contigs.sam")
    shell:
        "minimap2 --MD -ax sr {input.contigs} {input.filtered_reads} 1> {output.mapped} 2> {log.err_file}"


rule_name = "post_assembly__samtools_stats"
rule post_assembly__samtools_stats:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        mapped = rules.post_assembly__mapping.output.mapped
    output:
        stats = f"{component['name']}/contigs.stats",
    shell:
        "samtools stats {input.mapped} 1> {output.stats} 2> {log.err_file}"


rule_name = "post_assembly__pileup"
rule post_assembly__pileup:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        mapped = rules.post_assembly__mapping.output.mapped
    output:
        coverage = temp(f"{component['name']}/contigs.cov"),
        pileup = f"{component['name']}/contigs.pileup"
    shell:
        "pileup.sh in={input.mapped} basecov={output.coverage} out={output.pileup} 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__depth"
rule summarize__depth:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        coverage = rules.post_assembly__pileup.output.coverage
    params:
        component_json = component.json
    output:
        _file = f"{component['name']}/contigs.sum.cov",
        _file2= f"{component['name']}/contigs.bin.cov"
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "rule__summarize_depth.py")


rule_name = "post_assembly__call_variants"
rule post_assembly__call_variants:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__skesa.output,
        mapped = rules.post_assembly__mapping.output.mapped
    output:
        variants = temp(f"{component['name']}/contigs.vcf"),
    shell:
        "callvariants.sh in={input.mapped} vcf={output.variants} ref={input.contigs} ploidy=1 clearfilters 1> {log.out_file} 2> {log.err_file}"


rule_name = "summarize__variants"
rule summarize__variants:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        variants = rules.post_assembly__call_variants.output.variants
    params:
        component_json = component.json
    output:
        _file = f"{component['name']}/contigs.variants",
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "rule__summarize_variants.py")


rule_name = "rename_contigs"
rule rename_contigs:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    # Dynamic
    input:
        contigs = rules.assembly__skesa.output,
    output:
        contigs = f"{component['name']}/{sample['display_name']}.fasta"
    params:
        sample_name = sample['display_name']
    shell:
        "sed -e 's/Contig/{params.sample_name}/' {input.contigs} > {output.contigs}"
#* Dynamic section: end ****************************************************************************

#- Templated section: start ------------------------------------------------------------------------
rule_name = "datadump"
rule datadump:
    # Static
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        #* Dynamic section: start ******************************************************************
        rules.rename_contigs.output.contigs,  # Needs to be output of final rule
        rules.assembly_check__quast_on_contigs.output.quast,
        rules.post_assembly__samtools_stats.output.stats,
        rules.summarize__variants.output._file,
        rules.summarize__depth.output._file
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
