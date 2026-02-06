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

envvars:
    "BIFROST_INSTALL_DIR",
    "CONDA_PREFIX"

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

rule_name = "assembly__spades"
rule assembly__spades:
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
        rules.check_requirements.output.check_file,
        filtered_reads = sample['categories']['paired_reads']['summary']['trimmed']
    output:
        outputdir = Directory(f"{component['name']}/spades")
        scaffolds = f"{component['name']}/scaffolds.fasta"
    threads: 8
    shell:
        "spades -1 {input.filtered_reads[0]} -2 {input.filtered_reads[1]} --threads {threads} --isolate -o {output.outputdir} 1> {log.out_file} 2> {log.err_file} && cp {output.outputdir}/scaffolds.fasta {output.scaffolds}"


rule_name = "rename_scaffolds"
rule rename_scaffolds:
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
        scaffolds_in = rules.assembly__spades.output.scaffolds,
    output:
        scaffolds_out = f"{component['name']}/{sample['name']}.fasta"
    params:
        sample_name = sample['display_name']
    shell:
        "sed -e 's/NODE/{params.sample_name}/' {input.scaffolds_in} > {output.scaffolds_out}"

rule_name = "assembly_qc"
rule assembly_qc:
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
        scaffolds = rules.rename_scaffolds.output.scaffolds_out,
    output:
        scaffolds = f"{component['name']}/{sample['name']}_trimmed.fasta" 
        statistics = f"{component['name']}/{sample['name']}_trimmed_stat.tsv"
        failed_cov_scaffolds = f"{component['name']}/{sample['name']}_cov_fail.fasta"
        failed_cov_statistics = f"{component['name']}/{sample['name']}_cov_fail_stat.tsv"
        failed_length_scaffolds = f"{component['name']}/{sample['name']}_length_fail.fasta"
        failed_length_statistics = f"{component['name']}/{sample['name']}_length_fail_stat.tsv"
    params:
        sample_name = sample['display_name'],
        cov_threshold = 10
        min_length = 500
        passed_prefix=f"{component['name']}/{sample['name']}_trimmed"
        failed_cov_prefix=f"{component['name']}/{sample['name']}_cov_fail"
        failed_len_prefix=f"{component['name']}/{sample['name']}_length_fail"
    shell:
        """
        python rule__assembly_qc.py \
          --assembly "{input.scaffolds}" \
          --cov-threshold {params.cov_threshold} \
          --min-contig-len {params.min_length} \
          --passed {params.passed_prefix} \
          --failed-length {params.failed_len_prefix} \
          --failed-coverage {params.failed_cov_prefix} \
          --stdout \
          > "{log.out_file}" 2> "{log.err_file}"
        """

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
        #* Dynamic section: end ********************************************************************
    output:
        complete = rules.all.input
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")
#- Templated section: end --------------------------------------------------------------------------
