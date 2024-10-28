"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 10:30 AM 12/19/21
"""

import argparse
import logging
import os
import sys
import traceback
import time
from pprint import pprint

import decoil.output.metrics
import json
# import pickle5 as pickle
import pickle

# local modules
from decoil import __version__
import decoil.search.search as search
import decoil.encode.encode as encode
import decoil.encode.operations as operations
import decoil.search.quantify as quantify
import decoil.output.metrics as metrics
import decoil.search.cycles as cycles
import decoil.validate.compare as validate
from decoil.utils import QUAL
from decoil.utils import POS
from decoil.utils import PROG
from decoil.utils import VCF_PROP 

import decoil.output.parse as parse

log = logging.getLogger('reconstruct')
log.propagate = False
log.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)


def run_save_fragments(G, outfile):
	with open(outfile, "w") as f:
		f.write("#chr\tstart\tstop\tfid\tcoverage\n")
		fdict = G.get_fragments()
		for fid in fdict:
			_chr, start, end, coverage = fdict[fid].chr, fdict[fid].start, fdict[fid].end, fdict[fid].coverage
			if end == POS.MAX:
				if _chr not in ["chrM", "M"]:
					end = start + 2000
				else:
					end = start + 5
			f.write("""{}\t{}\t{}\t{}\t{}\n""".format(str(_chr), str(start), str(end), str(fid), str(int(coverage))))


def run_save_metrics(G, outputdir):
	"""
	Store metrics
	"""
	# get length and coverage fragment distribution
	glen = metrics.fragment_length_distribution(G)
	gcov = metrics.fragment_coverage_distribution(G)
	metrics.fragment_plot(glen, gcov, outputdir)

	# get mean coverage
	print("Mean coverage across sequenced regions", QUAL.MEAN_COVERAGE_WGS)


def run_save_simple_circles(graph):
	"""
	Store intermediate step (all simple circles) - debugging step
	"""
	bedfiledebug = "reconstruct.bed_debug"
	linksfiledebug = "reconstruct.links.txt_debug"

	log.info("7. Convert simple circles (all) to bed")
	metahash = graph.get_hash_canonical()
	simple_canonical = [metahash[m] for m in metahash]
	parse.convert_cycles2bed(simple_canonical, bedfiledebug, graph)

	log.info("8. Convert bed to links")
	parse.convert_bed2links(bedfiledebug, linksfiledebug)
 
 
def run_save_clusters(clusters, outputdir):
	"""
	Save the clusters
	"""
	clusterfile = os.path.join(outputdir,"clusters.json")
	with open(clusterfile, "w", encoding="utf-8") as json_file:
		json.dump(clusters, json_file, indent=4)
			

def run_save_v2_new(likelycandidates, G, ref_genome, outputdir, sample):

	os.makedirs(outputdir, exist_ok=True)

	fastafile = "reconstruct.fasta"
	fastafile_ecdna = "reconstruct.ecDNA.fasta"
	fastafile_ecdna_filtered = "reconstruct.ecDNA.filtered.fasta"
	jsonfile = "reconstruct.json"
	bedfile = "reconstruct.bed"
	bedfile_ecdna = "reconstruct.ecDNA.bed"
	bedfile_ecdna_filtered = "reconstruct.ecDNA.filtered.bed"

	linksfile = "reconstruct.links.txt"
	linksfile_ecdna = "reconstruct.links.ecDNA.txt"
	linksfile_ecdna_filtered = "reconstruct.links.ecDNA.filtered.txt"

	# htmlfile = "reconstruct.html"
	# htmlfile_ecdna = "reconstruct.ecDNA.html"
	# htmlfile_ecdna_filtered = "reconstruct.ecDNA.filtered.html"
	graphfile = "graph.gpickle"
	# graphplotfile = "graph.png"
	# treeplotfile = "tree.png"
	summaryfile = "summary.txt"
	# viz.plot_graph(G_tree, treeplotfile)

	log.info("6.1 Write to json")
	with open(jsonfile, "w", encoding='utf-8') as f:
		json.dump(likelycandidates, f, ensure_ascii=False, indent=4)

	log.info("6.2 Write to pickle format")
	with open(graphfile, 'wb') as f:
		pickle.dump(G, f, protocol=5)

	# log.info("6.3 Plot graph")
	# visualize.plot_graph(G, graphplotfile)

	log.info("6.3 Convert path to bed")
	parse.convert_path2bed_v2(likelycandidates, bedfile, keep=None)
	parse.convert_path2bed_v2(likelycandidates, bedfile_ecdna, keep="ecDNA", score_threshold=0)
	parse.convert_path2bed_v2(likelycandidates, bedfile_ecdna_filtered, keep="ecDNA", score_threshold=QUAL.FILTER_SCORE)

	log.info("6.4 Convert bed path to fasta")
	parse.convert_bed2fasta(bedfile, fastafile, ref_genome, version=2)
	parse.convert_bed2fasta(bedfile_ecdna, fastafile_ecdna, ref_genome, version=2)
	parse.convert_bed2fasta(bedfile_ecdna_filtered, fastafile_ecdna_filtered, ref_genome, version=2)

	log.info("6.5 Convert bed to links - required for visualization")
	parse.convert_bed2links(bedfile, linksfile)
	parse.convert_bed2links(bedfile_ecdna, linksfile_ecdna)
	parse.convert_bed2links(bedfile_ecdna_filtered, linksfile_ecdna_filtered)

	log.info("6.6. Summary of cycles with labels and annotation")
	parse.create_summary(likelycandidates, summaryfile)


def debug_cleaning(G):
	frags = G.get_fragments()
	for chr_ in ["chr1", "chr2", "chr3", "chr6"]:
		lenfrags = len([k for k in frags if frags[k].chr == chr_])
		lenintervals = len(G.get_fragment_intervals_bychr(chr_))
		print("After cleaning, remaining fragments vs fragments intervals", chr_, lenfrags, lenintervals)
	print("After cleaning fragments ", len(G.get_fragments()))
	print(G)


def run_reconstruction(vcffile, bigwigfile, bamfile, outputdir, ref_genome,
					   name="default_circle", svcaller=VCF_PROP.SNIFFLES1, fast=False):
	"""
	Reconstruct circle by finding the longest circular path in the graph
	"""
	os.makedirs(outputdir, exist_ok=True)
	# change to local directory
	currpath = os.getcwd()
	os.chdir(outputdir)

	# 0. Set thesholds
	metrics.set_wgs(bigwigfile)
	metrics.set_min_fragment_coverage(QUAL.MEAN_COVERAGE_WGS, QUAL.MINIMAL_FRAGMENT_COVERAGE)

	# 1 Graph generation
	# 1.1 Encode
	G = encode.run_encoding(vcffile, bigwigfile, bamfile, outputdir, svcaller=svcaller)
	run_save_fragments(G, "fragments_initial.bed")

	# 1.3 Add spatial
	G = operations.add_spatial_edges_new(G, bamfile)

	# 2. Clean
	# 2.1 Remove standalone fragments
	G = operations.remove_standalone_fragments(G)
	run_save_fragments(G, "fragments_clean_1.bed")

	# 2.2 Remove low coverage fragments
	# dynamically setup minimal fragment coverage (better not use this)
	# metrics.set_threshold()
	G = operations.remove_lowcoverage_fragments(G, threshold=QUAL.MINIMAL_FRAGMENT_COVERAGE)
	run_save_fragments(G, "fragments_clean_2.bed")

	# 2.3 Remove short fragments
	G = operations.remove_short_fragments(G, threshold=QUAL.MINIMAL_FRAGMENT_SIZE)
	run_save_fragments(G, "fragments_clean_3.bed")

	# 2.4 Remove duplicated edges
	G = operations.remove_duplicated_edges(G)

	# 2.3 Clean low coverage / small size fragments
	# G = operations.run_cleaning(G, outputdir)
	# debug_cleaning(G)

	# 3. Metrics
	metrics.compute_normalized_coverage(G)
	metrics.set_max_coverage(G)

	# 4. Search circular paths
	# 4.1 Search for all simple circular paths in the tree
	simple_cycles = search.run_searching(G, outputdir, name)
	run_save_simple_circles(G)

	# 4.2 Combine paths
	search.run_combine_paths(simple_cycles, G)

	# 5. Quantify likely true circles
	likely_candidates, clusters = quantify.run_quantification(G)
	likely_candidates_filtered = cycles.filter_reconstructions(likely_candidates, G)
	# likely_candidates_filtered_secondround = quantify.run_lasso_second(likely_candidates_filtered, G)

	likely_candidates_transformed = cycles.transform_fid2genomicpos(likely_candidates_filtered, G)
	likely_candidates_simplified = cycles.simplify_likely_candidates(likely_candidates_transformed)

	# 7. Label and annotate
	candidates_annotated = metrics.annotate(likely_candidates_simplified)

	# 8. Save
	# run_save_metrics(G, outputdir)
	# likelycandidates, G, ref_genome, outputdir, bigwigfile, gtffile, sample
	run_save_v2_new(candidates_annotated, G, ref_genome, outputdir, name)
	os.chdir(currpath)
 
 
def process_commandline_decoil_only(subparsers):
	"""
	Commandline for Decoil (without the abfront processing of the bam).
	Assumes the user is experienced.
	"""
	parser_a = subparsers.add_parser(PROG.RECONSTRUCT,
									 help='Reconstruct ecDNA',
									 usage='''decoil reconstruct -b <bamfile> -i <vcffile> -c <coveragefile> --outputdir <outputdir> --name <sample> -r <reference_genome>''')
	requiredNamed = parser_a.add_argument_group('required named arguments')
	requiredNamed.add_argument('-b', '--bam', help='Bam file', required=True)
	requiredNamed.add_argument('-c', '--coverage', help='Coverage file (bigwig)', required=True)
	requiredNamed.add_argument('-i', '--vcf', help='Vcf file', required=True)
	requiredNamed.add_argument('-o', '--outputdir', help='Output directory', required=True)
	requiredNamed.add_argument('--name', help='Name of the sample', required=True)
	requiredNamed.add_argument('-r', '--reference-genome', help='Reference genome (fasta)', required=True)
	parser_a.add_argument('-g', '--annotation-gtf', help='GTF annotation', required=False)
	parser_a.add_argument('-d', '--debug', help='Debug mode', action='count', default=0)
	parser_a.add_argument('--fast', help='Reconstruct fast (not accurate and does not require a bam file)',
						  action='count', default=0)

	# optional parameters
	parser_a.add_argument('--min-sv-len', help='Minimal SV length (default: %(default)sX)',
						  required=False, default=QUAL.MINIMAL_SV_LEN, type=int)
	parser_a.add_argument('--fragment-min-cov', help='Minimal fragment coverage (default: %(default)sX)',
						  required=False, default=QUAL.MINIMAL_FRAGMENT_COVERAGE, type=int)
	parser_a.add_argument('--fragment-min-size', help='Minimal fragment size (default: %(default)sbp)',
						  required=False, default=QUAL.MINIMAL_FRAGMENT_SIZE, type=int)
	parser_a.add_argument('--min-vaf', help='Minimal VAF acceptance SV (default: %(default)s)',
						  required=False, default=QUAL.MIN_VAF, type=float)
	parser_a.add_argument('--min-cov-alt', help='Minimal supporting reads SV (default: %(default)sX)',
						  required=False, default=QUAL.MIN_COV_ALT, type=int)
	parser_a.add_argument('--max-explog-threshold', help='Maximal score; better not change this (default: %(default)s)',
						  required=False, default=QUAL.EXPLOG_THRESHOLD, type=float)
	parser_a.add_argument('--min-cov', help='Minimal coverage on site (default: %(default)sX)',
						  required=False, default=QUAL.MIN_COV, type=int)
	parser_a.add_argument('--sv-caller', help="""SV caller name matching the VCF input {}""".format(VCF_PROP.SVCALLERS), required=False,
						  default=VCF_PROP.SNIFFLES1)
	parser_a.add_argument('--filter-score',
						  help='Filter circular structures by estimated copy-number (default: %(default)d)',
						  required=False, default=0, type=int)
	parser_a.set_defaults(which=PROG.RECONSTRUCT)

	parser_b = subparsers.add_parser(PROG.VALIDATE, help='Validate true against reconstruct ecDNA')
	parser_b.add_argument('-a', help='Configuration file for the true ecDNA', required=True)
	parser_b.add_argument('-b', help='Configuration file for the reconstructed ecDNA', required=True)
	parser_b.add_argument('-o', '--outputfile', help='Output file', required=True)
	parser_b.add_argument('-r', '--genome', help='Reference genome (fasta)', required=True)
	parser_b.add_argument('-t', '--temp', help='Temporary folder', required=False,
						  default='/data/gpfs-1/users/giurgium_c/scratch/temp')
	parser_b.set_defaults(which=PROG.VALIDATE)

	parser_c = subparsers.add_parser(PROG.FILTER, help='Filter VCF file')
	parser_c.add_argument('-i', '--vcf', help='Input vcf', required=True)
	parser_c.add_argument('-o', '--outputdir', help='Output directory', required=True)
	parser_c.set_defaults(which=PROG.FILTER)
 
	parser_d = subparsers.add_parser(PROG.CHECK, help='Check integrity of files')
	parser_d.add_argument('-i', '--input', help='Input file', required=True)
	parser_d.add_argument('-x', '--format', help='File format {VCF}', required=True, default='VCF')
	parser_d.set_defaults(which=PROG.CHECK)

	return subparsers


def process_commandline_decoil_fullpipeline(parser, subparsers):
	parser.add_argument('-n', '--dry-run', action='store_true')
	parser.add_argument('-f', '--force', action='store_true')
	parser.add_argument('-c', '--use-conda', action='store_true')

	# sv calling workflow
	parser_a = subparsers.add_parser(PROG.SV_ONLY,
									 description='Perform preprocessing (sv calling + coverage track)',
									 usage='''decoil-pipeline sv-only --bam <input> --outputdir <outputdir> --name <sample>''')
	requiredNamed = parser_a.add_argument_group('required named arguments')
	requiredNamed.add_argument('-b', '--bam', required=True, type=str)
	requiredNamed.add_argument('-o', '--outputdir', required=True, type=str)
	requiredNamed.add_argument('--name', help='Name of the sample', required=True)
	parser_a.add_argument('--threads', help='Number of threads (default: %(default)s)', required=False, default=4, type=int)
	parser_a.set_defaults(which=PROG.SV_ONLY)

	# sv-reconstruct workflow
	parser_b = subparsers.add_parser(PROG.SV_RECONSTRUCT,
									 description='Perform preprocessing (sv calling + coverage track) and reconstruction',
									 usage='''decoil-pipeline sv-reconstruct --bam <input> --outputdir <outputdir> --name <sample> -r <reference-genome>''')
	requiredNamed = parser_b.add_argument_group('required named arguments')
	requiredNamed.add_argument('-b', '--bam', required=True, type=str)
	requiredNamed.add_argument('-o', '--outputdir', required=True, type=str)
	requiredNamed.add_argument('--name', help='Name of the sample', required=True)
	requiredNamed.add_argument('-r', '--reference-genome', help='Reference genome (fasta)', required=True)
	parser_b.add_argument('--threads', help='Number of threads (default: %(default)s)', required=False, default=4, type=int)
	parser_b.add_argument('-g', '--annotation-gtf', help='GTF annotation', required=False)
	parser_b.add_argument('-d', '--debug', help='Debug mode', action='count', default=0)
	# optional parameters
	parser_b.add_argument('--min-sv-len', help='Minimal SV length (default: %(default)sX)',
						  required=False, default=QUAL.MINIMAL_SV_LEN, type=int)
	parser_b.add_argument('--fragment-min-cov', help='Minimal fragment coverage (default: %(default)sX)',
						  required=False, default=QUAL.MINIMAL_FRAGMENT_COVERAGE, type=int)
	parser_b.add_argument('--fragment-min-size', help='Minimal fragment size (default: %(default)sbp)',
						  required=False, default=QUAL.MINIMAL_FRAGMENT_SIZE, type=int)
	parser_b.add_argument('--min-vaf', help='Minimal VAF acceptance SV (default: %(default)s)',
						  required=False, default=QUAL.MIN_VAF, type=float)
	parser_b.add_argument('--min-cov-alt', help='Minimal supporting reads SV (default: %(default)sX)',
						  required=False, default=QUAL.MIN_COV_ALT, type=int)
	parser_b.add_argument('--max-explog-threshold',
						  help='Maximal score; better not change this (default: %(default)s)',
						  required=False, default=QUAL.EXPLOG_THRESHOLD, type=float)
	parser_b.add_argument('--min-cov', help='Minimal coverage on site (default: %(default)sX)',
						  required=False, default=QUAL.MIN_COV, type=int)
	parser_b.add_argument('--filter-score',
						  help='Filter circular structures by estimated proportions (default: %(default)d)',
						  required=False, default=0, type=int)
	parser_b.set_defaults(which=PROG.SV_RECONSTRUCT)

	# reconstruct-only workflow
	parser_d = subparsers.add_parser(PROG.RECONSTRUCT_ONLY,
									 description='Perform reconstruction only',
									 usage='''decoil-pipeline reconstruct-only --bam <input> --outputdir <outputdir> --name <sample> -r <reference-genome>''')
	requiredNamed = parser_d.add_argument_group('required named arguments')
	requiredNamed.add_argument('-b', '--bam', help='Bam file', required=True, type=str)
	requiredNamed.add_argument('-o', '--outputdir', required=True, type=str)
	requiredNamed.add_argument('--name', help='Name of the sample', required=True)
	requiredNamed.add_argument('-r', '--reference-genome', help='Reference genome (fasta)', required=True)
	parser_d.add_argument('--threads', help='Number of threads (default: %(default)s)', required=False, default=4, type=int)
	parser_d.add_argument('-g', '--annotation-gtf', help='GTF annotation', required=False)
	parser_d.add_argument('-d', '--debug', help='Debug mode', action='count', default=0)

	# optional parameters
	parser_d.add_argument('--min-sv-len', help='Minimal SV length (default: %(default)sX)',
						  required=False, default=QUAL.MINIMAL_SV_LEN, type=int)
	parser_d.add_argument('--fragment-min-cov', help='Minimal fragment coverage (default: %(default)sX)',
						  required=False, default=QUAL.MINIMAL_FRAGMENT_COVERAGE, type=int)
	parser_d.add_argument('--fragment-min-size', help='Minimal fragment size (default: %(default)sbp)',
						  required=False, default=QUAL.MINIMAL_FRAGMENT_SIZE, type=int)
	parser_d.add_argument('--min-vaf', help='Minimal VAF acceptance SV (default: %(default)s)',
						  required=False, default=QUAL.MIN_VAF, type=float)
	parser_d.add_argument('--min-cov-alt', help='Minimal supporting reads SV (default: %(default)sX)',
						  required=False, default=QUAL.MIN_COV_ALT, type=int)
	parser_d.add_argument('--max-explog-threshold',
						  help='Maximal score; better not change this (default: %(default)s)',
						  required=False, default=QUAL.EXPLOG_THRESHOLD, type=float)
	parser_d.add_argument('--min-cov', help='Minimal coverage on site (default: %(default)sX)',
						  required=False, default=QUAL.MIN_COV, type=int)
	parser_d.add_argument('--filter-score',
						  help='Filter circular structures by estimated proportions (default: %(default)d)',
						  required=False, default=0, type=int)
	parser_d.set_defaults(which=PROG.RECONSTRUCT_ONLY)

	return parser, subparsers


def process_commandline(sysargs, pipeline=False):
	
	if len(sysargs) == 0:
		sysargs.append("--help")
  
	parser = argparse.ArgumentParser(prog=PROG.DECOIL,
									 description="""Decoil {}: reconstruct ecDNA from long-read data""".format(
										 decoil.__version__),
									 usage='''
          decoil-pipeline [options] <run-mode> <parameters> [<target>]
          or
          decoil <run-mode> <paramters>''')
	parser.add_argument('--version', action='version',
						version='%(prog)s {}'.format(decoil.__version__))
	subparsers = parser.add_subparsers(help='sub-command help')

	# if decoil run as full pipeline
	if pipeline:
		parser, subparsers = process_commandline_decoil_fullpipeline(parser, subparsers)
	else:
		# decoil run as a standalone tool
		subparsers = process_commandline_decoil_only(subparsers)

	# read arguments
	args = parser.parse_args(sysargs)
	if args:
		subcommand = args.which

	return subcommand, args, parser


def setup_defaults(args):
	QUAL.MINIMAL_FRAGMENT_COVERAGE = args.fragment_min_cov
	QUAL.MINIMAL_FRAGMENT_SIZE = args.fragment_min_size
	QUAL.MINIMAL_SV_LEN = args.min_sv_len
	QUAL.MIN_VAF = args.min_vaf
	QUAL.MIN_COV = args.min_cov
	QUAL.MIN_COV_ALT = args.min_cov_alt
	QUAL.FILTER_SCORE = args.filter_score
	QUAL.EXPLOG_THRESHOLD = args.max_explog_threshold


def main(sysargs=sys.argv[1:]):

	try:
		start_time = time.time()

		if len(sysargs) == 0:
			sysargs.append('--help')
  
		# parse arguments
		subcommand, args, parser = process_commandline(sysargs)

		if subcommand == PROG.VALIDATE:
			# start validation
			validate.compare_true_reconstruct(os.path.abspath(args.a),
											  os.path.abspath(args.b),
											  os.path.abspath(args.reference_genome),
											  os.path.abspath(args.outputfile),
											  temp=os.path.abspath(args.temp))

		elif subcommand == PROG.RECONSTRUCT:
			# start reconstruction

			# start debug mode or not
			if args.debug == 1:
				log.setLevel(logging.DEBUG)
			else:
				log.setLevel(logging.INFO)

			# setup configuration parameters
			setup_defaults(args)

			# not use bam file
			fast = None
			bam = None
			if args.fast == 1:
				bam = None
				fast = True
			else:
				bam = os.path.abspath(args.bam)
				fast = False

			run_reconstruction(os.path.abspath(args.vcf),
							   os.path.abspath(args.coverage),
							   bam,
							   os.path.abspath(args.outputdir),
							   os.path.abspath(args.reference_genome),
							   name=args.name,
							   svcaller=args.sv_caller,
							   fast=fast)

		elif subcommand == PROG.FILTER:
			# start filtering
			operations.filter(os.path.abspath(args.vcf), os.path.abspath(args.outputdir))

		print("-------")
		print("Status: Successfully finished")
		print("User time (seconds): decoil", round(time.time() - start_time, 2))
		print("Subprogram: decoil", subcommand)
		print("Params:", args)
		print("#######")
		print()

	except AttributeError:
		parser.print_help()
		traceback.print_exc()
	except Exception:

		print("-------")
		print("Status: Failed")
		print("User time (seconds):", round(time.time() - start_time, 2))
		print("Subprogram: decoil", subcommand)
		print("Params:", args)
		print("#######")
		print()
		traceback.print_exc()


if __name__ == '__main__':
	main()

