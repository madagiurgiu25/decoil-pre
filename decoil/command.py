#! /usr/bin/env python
"""
Execution script for snakemake workflows.
https://github.com/ctb/2018-snakemake-cli/blob/master/run
"""

import argparse
import os
import os.path
import snakemake
import pprint
import json
import os
import sys
import traceback
import time
from pprint import pprint

# local modules
import decoil
from decoil import main as dp
from decoil.utils import QUAL
from decoil.utils import POS
from decoil.utils import PROG
from decoil.utils import VCF_PROP

# from . import _program

thisdir = os.path.abspath(os.path.dirname(__file__))
parentdir = os.path.join(thisdir, "")
cwd = os.getcwd()


def check_workflow(workflow):

    # first, find the Snakefile
    snakefile = os.path.join(thisdir, "cli", "Snakefile")
    if not os.path.exists(snakefile):
        sys.stderr.write("Error: cannot find Snakefile at {}\n".format(snakefile))
        sys.exit(-1)

    # next, find the workflow config file
    workflowfile = None
    if os.path.exists(workflow) and not os.path.isdir(workflow):
        workflowfile = workflow
    else:
        for suffix in ("", ".json"):
            tryfile = os.path.join(thisdir, "cli", workflow + suffix)
            if os.path.exists(tryfile) and not os.path.isdir(tryfile):
                sys.stderr.write("Found workflowfile at {}\n".format(tryfile))
                workflowfile = tryfile
                break

    if not workflowfile:
        sys.stderr.write("Error: cannot find workflowfile {}\n".format(workflow))
        sys.exit(-1)

    # find workflow path
    with open(workflowfile, "rt") as fp:
        workflow_info = json.load(fp)
    target = workflow_info["workflow_target"]
    configfile = workflow_info["workflow_config"]

    # read configuration
    with open(os.path.join(thisdir, "cli", configfile)) as f:
        config = json.load(f)
    config["workflow"] = workflow

    # load params if exists
    return snakefile, workflowfile, target, config


def add_defaults_toconfigs(args, config):
    config["params"] = ""
    config["params"] += " --fragment-min-cov " + str(args.fragment_min_cov)
    config["params"] += " --fragment-min-size " + str(args.fragment_min_size)
    config["params"] += " --fragment-max-cov " + str(args.fragment_max_cov)
    config["params"] += " --min-sv-len " + str(args.min_sv_len)
    config["params"] += " --min-vaf " + str(args.min_vaf)
    config["params"] += " --min-cov-alt " + str(args.min_cov_alt)
    config["params"] += " --min-cov " + str(args.min_cov)
    config["params"] += " --max-explog-threshold " + str(args.max_explog_threshold)
    config["params"] += " --filter-score " + str(args.filter_score)
    config["params"] += " --sv-caller " + str(args.sv_caller)
    if args.extend_allowed_chr != "":
        config["params"] += " --extend-allowed-chr " + args.extend_allowed_chr
    if args.fast == 1:
        config["params"] += " --fast "
    if args.skip == 1:
        config["params"] += " --skip "
    if args.multi == 1:
        config["params"] += " --multi "

def add_configs_sv_workflow(args, config):
    config["bam"] = args.bam
    config["outputdir"] = args.outputdir
    config["name"] = args.name
    config["svcaller"] = args.sv_caller


def add_configs_reconstruct_only_workflow(args, config):
    config["bam"] = args.bam
    config["outputdir"] = args.outputdir
    config["name"] = args.name
    config["reference_genome"] = args.reference_genome
    config["annotation_gtf"] = args.annotation_gtf
    config["svcaller"] = args.sv_caller 
    config = add_defaults_toconfigs(args, config)

def add_configs_sv_reconstruct_workflow(args, config):
    config["bam"] = args.bam
    config["outputdir"] = args.outputdir
    config["name"] = args.name
    config["reference_genome"] = args.reference_genome
    config["annotation_gtf"] = args.annotation_gtf
    config["svcaller"] = args.sv_caller  # for full pipeline only sniffles1 allowed
    config["filt_version"] = args.filt_version
    config = add_defaults_toconfigs(args, config)

# def add_configs_sv_reconstruct_plot_workflow(args,config):
#     add_configs_sv_workflow(args,config)
#     config['plot'] = True
#     config['window'] = args.plot_window
#     config['filterscore'] = args.plot_filter_score
#     config['filtertop'] = args.plot_top

# def add_configs_plot_workflow(args,config):
#     config['outputdir'] = args.outputdir
#     config['name'] = args.name
#     config['reference_genome'] = args.reference_genome
#     config['annotation_gtf'] = args.annotation_gtf
#     config['plot'] = True
#     config['window'] = args.plot_window
#     config['filterscore'] = args.plot_filter_score
#     config['filtertop'] = args.plot_top


def entry_point(sysargs=sys.argv[1:]):

    try:
        start_time = time.time()
        subcommand, args, parser = dp.process_commandline(sysargs, pipeline=True)

        # check if workflow exists
        snakefile, workflowfile, target, config = check_workflow(subcommand)

        config["container"] = False
        outputdir = os.path.abspath(args.outputdir)

        os.makedirs(outputdir, exist_ok=True)
        config["threads"] = args.threads

        if subcommand == PROG.SV_ONLY:
            add_configs_sv_workflow(args, config)

        if subcommand == PROG.RECONSTRUCT_ONLY:
            add_configs_reconstruct_only_workflow(args, config)

        if subcommand == PROG.SV_RECONSTRUCT:
            add_configs_sv_reconstruct_workflow(args, config)

        # if subcommand == PROG.SV_RECONSTRUCT_PLOT:
        #     add_configs_sv_reconstruct_plot_workflow(args, config)

        # if subcommand == PROG.PLOT_ONLY:
        #     add_configs_plot_workflow(args, config)

        # if subcommand in [PROG.SV_RECONSTRUCT, PROG.SV_RECONSTRUCT_PLOT, PROG.RECONSTRUCT_ONLY]:
        if subcommand in [PROG.SV_RECONSTRUCT, PROG.RECONSTRUCT_ONLY]:

            dp.setup_defaults(args)
            config["decoil"] = decoil.__version__

        print("--------")
        print("details!")
        print("\tsnakefile: {}".format(snakefile))
        print("\tworkflowfile: {}".format(workflowfile))
        print("\tparams: {}".format(config))
        print("\ttarget: {}".format(target))
        print("--------")

        # run workflows
        status = snakemake.snakemake(
            snakefile,
            targets=[target],
            printshellcmds=True,
            dryrun=args.dry_run,
            forceall=args.force,
            use_conda=False,
            config=config,
            cores=int(args.threads)
        )

        if status == False:
            raise Exception("Snakemake failed")
        print("-------")
        print("Status: Successfully finished")
        print("User time (seconds): decoil", round(time.time() - start_time, 2))
        print("#######")
        print()

    except AttributeError:
        traceback.print_exc()
        parser.print_help()
    except Exception:
        print("-------")
        print("Status: Failed")
        print("User time (seconds):", round(time.time() - start_time, 2))
        print("#######")
        print()
        traceback.print_exc()


if __name__ == "__main__":

    status = entry_point()

# python setup.py build install
#
