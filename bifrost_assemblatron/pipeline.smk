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
from snakemake.io import directory
import datetime

os.umask(0o2)

try:
    sample_ref = SampleReference(_id=config.get('sample_id', None), name=config.get('sample_name', None))
    sample: Sample = Sample.load(sample_ref)
    if sample is None:
        raise Exception("invalid sample passed")

    component_ref = ComponentReference(name=config['component_name'])
    component: Component = Component.load(reference=component_ref)
    if component is None:
        raise Exception("invalid component passed")

    samplecomponent_ref = SampleComponentReference(
        name=SampleComponentReference.name_generator(sample.to_reference(), component.to_reference())
    )
    samplecomponent = SampleComponent.load(samplecomponent_ref)
    if samplecomponent is None:
        samplecomponent = SampleComponent(
            sample_reference=sample.to_reference(),
            component_reference=component.to_reference()
        )

    common.set_status_and_save(sample, samplecomponent, "Running")

except Exception:
    print(traceback.format_exc(), file=sys.stderr)
    raise Exception("failed to set sample, component and/or samplecomponent")

onerror:
    if not samplecomponent.has_requirements():
        common.set_status_and_save(sample, samplecomponent, "Requirements not met")
    if samplecomponent["status"] == "Running":
        common.set_status_and_save(sample, samplecomponent, "Failure")

envvars:
    "BIFROST_INSTALL_DIR",
    "CONDA_PREFIX"

# -------------------------------------------------------------------------
# MAIN + TIMING
# -------------------------------------------------------------------------

rule all:
    input:
        f"{component['name']}/datadump_complete"
    run:
        common.set_status_and_save(sample, samplecomponent, "Success")

rule set_time_start:
    output:
        start_file = temp(f"{component['name']}/time_start.txt")
    run:
        import time
        with open(output.start_file, "w") as fh:
            fh.write(str(time.time()))

rule setup:
    input:
        rules.set_time_start.output.start_file
    output:
        init_file = touch(temp(f"{component['name']}/initialized"))
    run:
        samplecomponent["path"] = os.path.join(os.getcwd(), component["name"])
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
        folder = rules.setup.output.init_file
    output:
        check_file = f"{component['name']}/requirements_met"
    params:
        samplecomponent
    run:
        if samplecomponent.has_requirements():
            with open(output.check_file, "w") as fh:
                fh.write("")

#* Dynamic section: start **************************************************************************

rule_name = "assembly__spades"
rule assembly__spades:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.check_requirements.output.check_file,
        filtered_reads = sample["categories"]["trimmed_reads"]["summary"]["data"]
    output:
        outputdir = directory(f"{component['name']}/spades"),
        scaffolds = f"{component['name']}/scaffolds.fasta",
        threads_file = temp(f"{component['name']}/threads_used.txt"),
        tool_version = temp(f"{component['name']}/tool_version.txt")
    params:
        threads = 10
    shell: r"""
        spades.py -1 {input.filtered_reads[0]} -2 {input.filtered_reads[1]} -t {params.threads} --isolate -o {output.outputdir} 1> {log.out_file} 2> {log.err_file}

        cp {output.outputdir}/scaffolds.fasta {output.scaffolds}

        echo {params.threads} > {output.threads_file}

        spades.py -v > {output.tool_version}	
        """

rule_name = "rename_scaffolds"
rule rename_scaffolds:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        scaffolds_in = rules.assembly__spades.output.scaffolds
    output:
        scaffolds_out = f"{component['name']}/{sample['name']}.fasta"
    params:
        sample_name = sample["display_name"]
    shell:
        "sed -e \"s/NODE/{params.sample_name}/\" {input.scaffolds_in} > {output.scaffolds_out}"

rule_name = "assembly_qc"
rule assembly_qc:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log",
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        scaffolds = rules.rename_scaffolds.output.scaffolds_out
    output:
        scaffolds = f"{component['name']}/{sample['name']}_trimmed.fasta",
        statistics = f"{component['name']}/{sample['name']}_trimmed_stat.tsv",
        failed_cov_scaffolds = f"{component['name']}/{sample['name']}_cov_fail.fasta",
        failed_cov_statistics = f"{component['name']}/{sample['name']}_cov_fail_stat.tsv",
        failed_length_scaffolds = f"{component['name']}/{sample['name']}_length_fail.fasta",
        failed_length_statistics = f"{component['name']}/{sample['name']}_length_fail_stat.tsv",
    params:
        sample_name = sample["display_name"],
        cov_threshold = 10,
        min_length = 500,
        passed_prefix = f"{component['name']}/{sample['name']}_trimmed",
        failed_cov_prefix = f"{component['name']}/{sample['name']}_cov_fail",
        failed_len_prefix = f"{component['name']}/{sample['name']}_length_fail",
        qc_script = os.path.join(os.path.dirname(workflow.snakefile), "rule__assembly_qc.py")
    shell:
        r"""
        python {params.qc_script} \
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

# -------------------------------------------------------------------------
# END TIME + RUNTIME (FILE-BASED)
# -------------------------------------------------------------------------

rule set_time_end:
    input:
        rules.assembly_qc.output.scaffolds
    output:
        end_file = temp(f"{component['name']}/time_end.txt")
    run:
        import time
        with open(output.end_file, "w") as fh:
            fh.write(str(time.time()))

rule dump_info:
    input:
        start_file = rules.set_time_start.output.start_file,
        end_file = rules.set_time_end.output.end_file,
        threads_file = rules.assembly__spades.output.threads_file,
        spades_version = rules.assembly__spades.output.tool_version
    output:
        runtime_flag = temp(f"{component['name']}/runtime_set")
    run:
        import time
        from bifrostlib.datahandling import SampleComponent

        with open(input.start_file) as fh:
            t_start = float(fh.read().strip())
        with open(input.end_file) as fh:
            t_end = float(fh.read().strip())
        with open(input.threads_file) as fh:
            threads_used = int(fh.read().strip())
        with open(input.spades_version) as fh:
            spades_version = str(fh.read())

        runtime_minutes = (t_end - t_start) / 60.0
        print(f"runtime in minutes {runtime_minutes}")

        sc = SampleComponent.load(samplecomponent.to_reference())
        sc["time_start"] = datetime.datetime.fromtimestamp(t_start).strftime("%Y-%m-%d %H:%M:%S")
        sc["time_end"] = datetime.datetime.fromtimestamp(t_end).strftime("%Y-%m-%d %H:%M:%S")
        sc["time_running"] = round(runtime_minutes, 3)
        sc["threads_used"] = threads_used
        sc["tool_version"] = spades_version

        sc.save()

        with open(output.runtime_flag, "w") as fh:
            fh.write("done")


# -------------------------------------------------------------------------
# DATADUMP
# -------------------------------------------------------------------------

rule_name = "datadump"
rule datadump:
    message:
        f"Running step:{rule_name}"
    log:
        out_file = f"{component['name']}/log/{rule_name}.out.log",
        err_file = f"{component['name']}/log/{rule_name}.err.log"
    benchmark:
        f"{component['name']}/benchmarks/{rule_name}.benchmark"
    input:
        rules.assembly_qc.output.scaffolds,
        rules.dump_info.output.runtime_flag
    output:
        complete = f"{component['name']}/datadump_complete"
    params:
        samplecomponent_ref_json = samplecomponent.to_reference().json
    script:
        os.path.join(os.path.dirname(workflow.snakefile), "datadump.py")

