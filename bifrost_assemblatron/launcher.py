#!/usr/bin/env python3
"""
Launcher file for accessing dockerfile commands
"""
import argparse
import json
import subprocess
import os
import sys
import traceback
from bifrostlib import datahandling

COMPONENT: dict = datahandling.load_yaml(os.path.join(os.path.dirname(__file__), 'config.yaml'))

def parser(args):
    """
    Arg parsing via argparse
    """
    description: str = (
        f"-Description------------------------------------\n"
        f"{COMPONENT['details']['description']}"
        f"------------------------------------------------\n\n"
        f"*Run command************************************\n"
        f"docker run \ \n"
        f" -e BIFROST_DB_KEY=mongodb://<user>:<password>@<server>:<port>/<db_name> \ \n"
        f" {COMPONENT['install']['dockerfile']} \ \n"
        f"************************************************\n"
    )
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--install',
                        action='store_true',
                        help='Install/Force reinstall component')
    parser.add_argument('-info', '--info',
                        action='store_true',
                        help='Provides basic information on component')
    parser.add_argument('-out', '--outdir',
                        default=".",
                        help='Output directory')
    parser.add_argument('-id', '--sample_id',
                        action='store',
                        type=str,
                        help='Sample ID of sample in bifrost, sample has already been added to the bifrost DB')

    try:
        options: argparse.Namespace = parser.parse_args(args)
    except:
        parser.print_help()
        sys.exit(0)

    return options

def run_program(args: argparse.Namespace):
    if not datahandling.check_db_connection_exists():
        message: str = (
            f"ERROR: Connection to DB not establised.\n"
            f"please ensure env variable BIFROST_DB_KEY is set and set properly\n"
        )
        print(message)
    else:
        print(datahandling.get_connection_info())

    if args.info:
        show_info()
    elif args.install:
        install_component()
    elif args.sample_id is not None:
        run_sample(args)


def show_info():
    """
    Shows information about the component
    """
    message: str = (
        f"Component: {COMPONENT['name']}\n"
        f"Version: {COMPONENT['version']}\n"
        f"Details: {json.dumps(COMPONENT['details'], indent=4)}\n"
        f"Requirements: {json.dumps(COMPONENT['requirements'], indent=4)}\n"
        f"Output files: {json.dumps(COMPONENT['db_values_changes']['files'], indent=4)}\n"
    )
    print(message)


def install_component():
    component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']])
    # if len(component) == 1:
    #     print(f"Component has already been installed")
    if len(component) > 1:
        print(f"Component exists multiple times in DB, please contact an admin to fix this in order to proceed")
    else:
        #HACK: Installs based on your current directory currently. Should be changed to the directory your docker/singularity file is
        #HACK: Removed install check so you can reinstall the component. Should do this in a nicer way
        COMPONENT['install']['path'] = os.path.os.getcwd()
        datahandling.post_component(COMPONENT)
        component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']])
        if len(component) != 1:
            print(f"Error with installation of {COMPONENT['name']} {len(component)}\n")
            exit()


def run_sample(args: object):
    """
    Runs sample ID through snakemake pipeline
    """
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    os.chdir(args.outdir)

    sample: list[dict] = datahandling.get_samples(sample_ids=[args.sample_id])
    component: list[dict] = datahandling.get_components(component_names=[COMPONENT['name']])
    if len(component) == 0:
        print(f"component not found in DB, would you like to install it (Y/N)?:")
        install_component()

    if len(sample) == 0:
        # Invalid sample id
        print(f"sample_id not found in DB")
        pass
    elif len(sample) != 1:
        print(f"Error with sample_id")
    elif len(component) != 1:
        print(f"Error with component_id")
    else:
        print(f"snakemake -s /bifrost_{COMPONENT['display_name']}/bifrost_{COMPONENT['display_name']}/pipeline.smk --config sample_id={str(sample[0]['_id'])} component_id={str(component[0]['_id'])}")
        try:
            process: subprocess.Popen = subprocess.Popen(
                f"snakemake -s /bifrost_{COMPONENT['display_name']}/bifrost_{COMPONENT['display_name']}/pipeline.smk --config sample_id={str(sample[0]['_id'])} component_id={str(component[0]['_id'])}",
                stdout=sys.stdout,
                stderr=sys.stderr,
                shell=True
            )
            process.communicate()
        except:
            print(traceback.format_exc())

def run():
    args: argparse.Namespace = parser(sys.argv[1:])
    run_program(args)

if __name__ == '__main__':
    run()
