"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 10:22 AM 11/14/22

Collect different metrics
"""
import logging
import os
import sys
import math

import matplotlib.pyplot as plt
import numpy as np
import pyBigWig
import seaborn as sns

from decoil.utils import QUAL
from decoil.utils import META_CONFORMATION as mc
from decoil.utils import TOPOLOGY as tp

handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def fragment_length_distribution(graph):
	"""
	Compute fragment length distribution for the graph
	"""
	frags = graph.get_fragments()
	frag_length = []
	for fid in frags:
		frag_length.append(int(abs(frags[fid].start - frags[fid].end)))
	print(frag_length)
	return frag_length


def fragment_coverage_distribution(graph):
	"""
	Compute fragment coverage distribution for the graph
	"""
	frags = graph.get_fragments()
	frag_cov = []
	for fid in frags:
		frag_cov.append(int(frags[fid].coverage))
	print(frag_cov)
	return frag_cov


def fragment_plot(glen, gcov, outputdir):
	"""
	Plot frag len/cov distribution
	"""

	# add pseudocount
	glen = np.array(glen) + 0.01
	gcov = np.array(gcov) + 0.01

	os.makedirs(outputdir, exist_ok=True)
	plot_glen = os.path.join(outputdir, "metrics_frag_len_distribution.png")
	plot_gcov = os.path.join(outputdir, "metrics_frag_cov_distribution.png")
	plot_sct = os.path.join(outputdir, "metrics_frag_len_cov_correlation.png")
	
	g = sns.histplot(glen, bins=50, kde=True, log_scale=True)
	plt.title('Fragment length distribution')
	g.figure.savefig(plot_glen, dpi=500)
	plt.close(g.figure)
	
	g1 = sns.histplot(gcov, bins=50, kde=True, log_scale=True)
	plt.title('Fragment coverage distribution')
	g1.figure.savefig(plot_gcov, dpi=500)
	plt.close(g1.figure)
	
	g2 = sns.scatterplot(x=glen, y=gcov)
	plt.title('Fragment length vs. coverage(weight)')
	g2.figure.savefig(plot_sct, dpi=500)
	plt.close(g2.figure)


def compute_meancov(bwfile):
	bw = pyBigWig.open(bwfile)
	
	mean_perchr = []
	len_perchr = []
	
	if bw.isBigWig():
		all_chroms = bw.chroms()
		for _chr in all_chroms:
			# select only primary chromosomes
			if (_chr.startswith("chr") or (not _chr.startswith("K") and not _chr.startswith("J"))) and \
					(not _chr.endswith("fix") or not _chr.endswith("alt")):
				mean_perchr.append(bw.stats(_chr, 0, all_chroms[_chr], type="mean")[0])
				len_perchr.append(all_chroms[_chr])
	else:
		raise ValueError("""{} needs to be a bigWig file""".format(bwfile))
	
	len_perchr = np.array(len_perchr) / sum(len_perchr)
	if len(mean_perchr) == 0:
		return 0
	return sum(np.array(mean_perchr) * np.array(len_perchr))


def set_wgs(bigwigfile):
	"""
	Compute WGS mean coverage
	"""
	QUAL.MEAN_COVERAGE_WGS = compute_meancov(bigwigfile)
	print("Set QUAL.MEAN_COVERAGE_WGS to ", QUAL.MEAN_COVERAGE_WGS)


def set_max_coverage(graph):
	"""
	Set the max coverage == max coverage across fragments
	"""
	max = 0
	fragments = graph.get_fragments()
	for f in fragments:
		if fragments[f].coverage > max:
			max = fragments[f].coverage

	QUAL.MAX_COVERAGE = max
	print("Set QUAL.MAX_COVERAGE to ", QUAL.MAX_COVERAGE)

def set_min_fragment_coverage(wgs, preset_min_coverage):
	"""
	Maxim between 2xWGS and the preset minimal fragment coverage (user defined)
	"""
	QUAL.MINIMAL_FRAGMENT_COVERAGE = max(QUAL.MEAN_COVERAGE_WGS,QUAL.MINIMAL_FRAGMENT_COVERAGE)
	print("Set QUAL.MINIMAL_FRAGMENT_COVERAGE to ", QUAL.MINIMAL_FRAGMENT_COVERAGE)

def set_threshold():
	"""
	Set threshold for keeping fragments in the graph
	"""
	QUAL.MINIMAL_FRAGMENT_COVERAGE = max(QUAL.MINIMAL_FRAGMENT_COVERAGE, QUAL.MEAN_COVERAGE_WGS * 4)
	print("Minimal coverage across fragments", QUAL.MINIMAL_FRAGMENT_COVERAGE)


def compute_normalized_coverage(graph):
	"""
	Compute normalized coverage for all fragments in the graph
	"""
	# print("Compute normalized coverage")
	for fid in graph.fragments:
		graph.fragments[fid].set_norm_cov(QUAL.MEAN_COVERAGE_WGS)
		# print(fid, graph.fragments[fid].get_norm_cov())


def annotate_size(cycle):
	"""
	Total size
	"""
	total_size = 0
	for (chr, start, stop, strand, frag, circid, cov) in cycle:
		total_size += abs(int(start)-int(stop))
	return total_size / 1000000

def annotate_chrs(cycle):
	"""
	All chr of origin
	"""
	lchr = []
	for (chr, start, stop, strand, frag, circid, cov) in cycle:
		lchr.append(chr)
	return ",".join(list(set(lchr)))

def annotate_label(size):
	"""
	Label simple by size
	"""

	if size > QUAL.ECDNA_MINSIZE:
		return "ecDNA"
	else:
		return ""

def get_topology(topology_details):
	"""
	Return topology based on the json details
	"""
	topology = 0
	# determine topology
	if topology_details[mc.N_FRAG_STR] == 1:
		topology = tp.SIMPLE_EXCISION
	else:
		if topology_details[mc.SMALL_DEL_STR] == 1:
			topology = tp.SIMPLE_EVENTS
		if topology_details[mc.INV_STR] == 1:
			topology = tp.SIMPLE_EVENTS
		if topology_details[mc.SMALL_DEL_STR] == 1 and topology_details[mc.INV_STR] == 1:
			topology = tp.MIXED_SIMPLE_EVENTS
		if topology_details[mc.MULTI_REGION_STR] == 1:
			topology = tp.MULTI_REGION_INTRA_CHR
		if topology_details[mc.INTERCHR_STR] == 1:
			topology = tp.MULTI_REGION_INTER_CHR
		if topology_details[mc.DUP_STR] == 1:
			topology = tp.SIMPLE_DUPLICATIONS
		if topology_details[mc.FOLDBACK_STR] == 1:
			topology = tp.FOLDBACKS
	return topology

def annotate_topology(cycle):
	"""
	Annotate the topology of the cycle.
	We define the following topologies:


	Metaconformation:
		N_FRAG - number of fragments
		SMALL_DEL - allow small deletions on the right and left side of the fragment
		DUP - allow simple duplication
		INV - allow inversions
		INTERCHR - allow fragments to originate from multiple chromosomes
		MULTI_REGION - allow fragments to originate from multiple regions on same chromosome
		FOLDBACK - allow foldbacks, which are here defined as two overlapping genomic fragments, which are immediately chained in the ecDNA template

	Topologies (based on the metainformation):

		Simple excisions:                1     0  0  0   0   0  0
		Simple events:                 [2-10] 0|1 0  0   0   0  0
									   [2-10]  0  0  1   0   0  0
		Mixed simple events:           [2-10]  1  0  1   0   0  0
		Multi-region intra-chromosomal:[2-10] 0|1 0 0|1  0   1  0
		Multi-region inter-chromosomal:[2-10] 0|1 0 0|1  1   1  0
		Simple duplications:           [2-10] 0|1 1 0|1 0|1 0|1 0
		Foldbacks:

	Arguments:
		cycle (list): [(chr, start, stop, strand, frag, circid, cov),]
	"""
	topology = None
	topology_details = {mc.N_FRAG_STR:0,
	                    # mc.LEN_STR:0,
	                    mc.SMALL_DEL_STR:0,
	                    mc.DUP_STR:0,
	                    mc.INV_STR:0,
	                    mc.INTERCHR_STR:0,
	                    mc.MULTI_REGION_STR:0,
	                    mc.FOLDBACK_STR:0
	                    }
	previous_fragment = None
	visited_fragments = []
	for (chr, start, stop, strand, frag, circid, cov) in cycle:

		cstart = min(int(stop),int(start))
		cstop = max(int(stop),int(start))
		clen = abs(cstart-cstop)

		# track nr fragments
		topology_details[mc.N_FRAG_STR] += 1

		# check if duplication
		if frag in visited_fragments:
			topology_details[mc.DUP_STR] = 1

		if previous_fragment:
			
			pstart = min(int(previous_fragment[1]), int(previous_fragment[2]))
			pstop = max(int(previous_fragment[1]), int(previous_fragment[2]))
			pstrand = previous_fragment[3]
			# print("previous", pstart, pstop, pstrand)
			# print("current", cstart, cstop, strand)
			plen = abs(pstart-pstop)
			# 10% of the previous fragment size

			del_threshold = min(5000, int(0.1 * plen + 0.1 * clen))
            #        del_threshold = 1000000
			# print("plen", plen)
			# print("del tresh", del_threshold)
			# print("diff",abs(cstart - pstop), abs(pstart - cstop))

			
			if chr == previous_fragment[0]:
				# check small deletions or multi-region or duplication
				if 0 <= abs(cstart - pstart) <= 50 and 0 <= abs(cstop - pstop) <= 50:
					topology_details[mc.DUP_STR] = 1
				elif 10 < (cstart - pstop) < del_threshold or 10 < (pstart - cstop) < del_threshold:
					topology_details[mc.SMALL_DEL_STR] = 1
				elif (cstart - pstop) >= del_threshold or (pstart - cstop) >= del_threshold:
					topology_details[mc.MULTI_REGION_STR] = 1

				# check inversions
				if strand != pstrand:
					topology_details[mc.INV_STR] = 1

				# check foldback
				if pstart < start < pstop or pstart < stop < pstop and not (0 <= abs(cstart - pstart) <= 50 and 0 <= abs(cstart - pstart) <= 50):
					topology_details[mc.FOLDBACK_STR] = 1
			else:
				# chimeric / multiple chr
				topology_details[mc.INTERCHR_STR] = 1

		previous_fragment = (chr, start, stop, strand, frag, circid, cov)
		visited_fragments.append(frag)
		
		# print(topology_details)
		# print()

	topology = get_topology(topology_details)

	return topology, topology_details

def annotate(candidates):
	"""
	Annotate circles:
	- size
	- complexity score
	- class
	- genes
	- label as ecDNA
	"""
	for cycle in candidates:
		print(candidates[cycle])
		candidates[cycle]["size"] = annotate_size(candidates[cycle]["conf"])
		candidates[cycle]["chrs"] = annotate_chrs(candidates[cycle]["conf"])
		candidates[cycle]["label"] = annotate_label(candidates[cycle]["size"])
		candidates[cycle]["topology"], candidates[cycle]["topology_detailed"] = annotate_topology(candidates[cycle]["conf"])
	return candidates
