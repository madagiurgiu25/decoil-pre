"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 5:30 PM 9/19/22

Different operations:
- VCF file cleaning:
	- Check quality of the variants
	- Polish breakpoints
	- Merge fuzzy breakpoints
- Graph cleaning:
	- Remove low coverage fragments
	- Remove singleton fragments
- Graph properties editing
	- Compute mean coverage per fragment
	- Compute edge weight connecting 2 neighboring fragments
"""
import logging
import math
import os
import numpy as np
from collections import defaultdict
from copy import deepcopy

import pyBigWig
import pysam
import vcfpy

import decoil.utils as utils
from decoil.encode.graph import MultiGraph
import decoil.encode.encode as encode
from decoil.output.metrics import handler, formatter
from decoil.utils import GRAPH_EDGE_PROP as gep
from decoil.utils import GRAPH_PROP as gp
from decoil.utils import QUAL
from decoil.utils import VCF_PROP as vp

log = logging.getLogger('reconstruct.clean')
log.propagate = False
log.setLevel(logging.DEBUG)

handler.setLevel(logging.DEBUG)
handler.setFormatter(formatter)
log.addHandler(handler)


def in_window(array, pos1, pos2, distance):
	"""
	Check if two breakpoints lie in close proximity
	"""
	if abs(array[pos1] - array[pos2]) <= distance:
		return True
	return False


def merge_near_breakpoints(collection_breakpoints, svinfo, distance=QUAL.DISTANCE):
	"""
	Merge brekpoints in the close proximity
	"""
	
	# for each chromosome
	for chr in collection_breakpoints:
		
		breakpoints_todelete = {}  # ref_breakpoint : [list of breakpoints in close proximity]
		breakpoints_partners = {}  # breakpoint close proximity : ref_breakpoint
		nrbreakpoints = len(collection_breakpoints[chr])
		
		# if you have more than 2 breakpoints per chromosome
		if nrbreakpoints >= 2:
			
			# merge if positions close together
			i = 0
			while i < nrbreakpoints - 1:
				reference_breakpoint = collection_breakpoints[chr][i]
				j = i + 1
				while j < nrbreakpoints:
					# next breakpoint is not in close proximity
					if in_window(collection_breakpoints[chr], i, j, distance):
						# breakpoint is in close proximity
						if reference_breakpoint not in breakpoints_todelete:
							breakpoints_todelete[reference_breakpoint] = []
						close_breakpoint = collection_breakpoints[chr][j]
						breakpoints_todelete[reference_breakpoint].append(close_breakpoint)
						breakpoints_partners[close_breakpoint] = reference_breakpoint
						log.info("Remove {}".format(close_breakpoint))
						
						j = j + 1
					else:
						# breakpoints are not in close proximity
						break
				
				# jump over all brekpoints in close proximity
				i = j
		
		for refbr in breakpoints_todelete:
			refkey = '{}{}{}'.format(chr, utils.SEPARATOR, str(refbr))
			
			# check if refkey in svinfo
			if refkey not in svinfo:
				svinfo[refkey] = []
			
			# remove all fuzzy breakpoints and transfer sv mate to the reference breakpoint
			for br in breakpoints_todelete[refbr]:
				# remove breakpoint from close proximity
				collection_breakpoints[chr].remove(br)
				collapsed_key = '{}{}{}'.format(chr, utils.SEPARATOR, str(br))
				
				# move all the sv attached to this breakpoint to the refkey
				if collapsed_key in svinfo:
					# transfer children
					for elem in svinfo[collapsed_key]:
						svinfo[refkey].append(elem)
					del svinfo[collapsed_key]
		
		# check if there are fuzzy sv mate breakpoints and convert them to reference breakpoint
		copy_svinfo = defaultdict(list)
		for key in svinfo:
			copy_svinfo[key] = []
			for index, sv_end in enumerate(svinfo[key]):
				if int(sv_end[2]) in breakpoints_partners and chr == sv_end[1]:
					# change chrom position
					new_sv = (
						sv_end[0],
						sv_end[1],
						str(breakpoints_partners[int(sv_end[2])]),
						sv_end[3],
						sv_end[4],
						sv_end[5],
						sv_end[6],
						sv_end[7])
					copy_svinfo[key].append(new_sv)
				else:
					copy_svinfo[key].append(sv_end)
		svinfo = deepcopy(copy_svinfo)
	
	for key in collection_breakpoints:
		collection_breakpoints[key] = sorted(list(set(collection_breakpoints[key])))
	
	return collection_breakpoints, svinfo


def pass_filter(record):
	"""
	Check if SV record pass the quality filter
	"""
	
	val = record.calls[0].data.get('DV')
	if val == None:
		raise Exception("DV has not value. Your VCF might not be genotyped. Rerun SV calling using --genotype")
	v = int(val)
	
	val = record.calls[0].data.get('DR')
	if val == None:
		raise Exception("DR has not value. Your VCF might not be genotyped. Rerun SV calling using --genotype")
	dr = int(val)
 
	cov = (dr + v)
	vaf = (v / (v + dr))

	# precise breakpoint
	# if vp.PRECISE not in record.INFO or not record.INFO[vp.PRECISE]:
	# 	return False

	if record.INFO[vp.SVTYPE] not in vp.SV_COLLECTION_STRICT:
		return False

	# filter low cov (< 10X)
	if cov < QUAL.MIN_COV:
		return False
	
	# filter low VAF
	if vaf < QUAL.MIN_VAF:
		return False
	
	# filter variants with < 4X
	if v < QUAL.MIN_COV_ALT:
		return False
	
	if 'PASS' not in record.FILTER and 'STRANDBIAS' not in record.FILTER:
		return False
	
	if record.INFO[vp.SVTYPE] not in vp.SV_COLLECTION:
		log.warning("SV unknown")
		log.warning(record)
		return False

	if str(record.INFO[vp.CHR2]) not in vp.ALLOWED_CHR or str(record.CHROM) not in vp.ALLOWED_CHR:
		return False

	if record.INFO[vp.SVTYPE] != vp.BND and abs(int(record.INFO[vp.SVLEN])) < QUAL.MINIMAL_SV_LEN:
		return False

	# filter based on a negative exponential curve
	if np.power(3,-cov*vaf/np.log(cov+0.001)) >= QUAL.EXPLOG_THRESHOLD:
		return False
	
	return True


def filter(vcf, vcfout, svcaller=vp.SNIFFLES1):
	"""
	Read vcf file and store all breakpoints.

	Return:
		(breakpoints, sv info)
	"""
        
	encode.cleanvcf(vcf,vcf + "_clean.vcf")
	reader = vcfpy.Reader.from_path(vcf + "_clean.vcf")
	# outfile = vcf.replace(".vcf", ".filtered.vcf")
	writer = vcfpy.Writer.from_path(vcfout, reader.header)
	
	for record in reader:
		try:
			# check quality filter
			if not pass_filter(record):
				# skip record as this is good
				continue
			writer.write_record(record)
		except Exception as e:
			print("WARNING incorrect", record)
			continue

def remove_duplicated_edges(graph):
	"""
	The graph is a multigraph, so between 2 nodes multiple edges are allowed (e.g. duplication of single fragment).
	However, there are cases where having multiple edges connecting the same 2 nodes does not make sense.
	In the case of small deletions (~50-200bp) the graph will include for the same 2 nodes one DEL and one Spatial edge.
	Task: remove these
	"""
	to_remove = []
	visited_edges = {}
	nodes = graph.get_nodes()
	edges = graph.get_edges()

	for u in nodes:
		for v in graph.graph[u]:
			# nodes do not belong to same fragment and have more than 2 connected edges
			if u != v and nodes[u].parent_fragment != nodes[v].parent_fragment and len(graph.graph[u][v]) > 1:

				w_edges = np.array([edges[e].weight for e in graph.graph[u][v]])
				idx_max = np.where(w_edges == w_edges.max())[0][0]

				# remove all other edges which dont high supporting reads
				for i, e in enumerate(graph.graph[u][v]):
					if e not in visited_edges and edges[e].svtype != gp.FRAGMENT:
						if i != idx_max:
							to_remove.append((e,u,v))
						visited_edges[e] = None

	for e,u,v in to_remove:
		print("remove edge ", e, " connecting nodes u,v ", u, v)
		graph.remove_edge(e,u,v)

	return graph

def remove_lowcoverage_fragments(graph, threshold=QUAL.MINIMAL_FRAGMENT_COVERAGE):
	"""
	Remove fragments which are covered < 4X (QUAL.MINIMAL_FRAGMENT_COVERAGE) and remove subgraph

	Args:
		graph (decoil.graph.Mutigraph): SV graph

	Returns:
		filtered graph including only fragments with coverage >= 4X
	"""
	
	# remove low coverage fragments
	fragments = graph.get_fragments()
	toremove = []
	for fid in fragments:
		frag = fragments[fid]
		if frag.coverage < threshold:
			toremove.append(fid)
	
	for f in toremove:
		# graph.rewire_fragment_neighbors(f)
		graph.remove_fragment(f)

	return graph


def remove_standalone_fragments(graph):
	"""
	Remove fragments which are not connected to any SV
	"""
	# remove low coverage fragments
	fragments = graph.get_fragments()
	toremove = []
	for fid in fragments:
		frag_obj = fragments[fid]
		u = frag_obj.tail
		v = frag_obj.head
		# check if one node is a loose end
		if graph.degree(u) == 0 and graph.degree(v) == 0:
			toremove.append(fid)
	
	for fid in toremove:
		graph.remove_fragment(fid)
	
	return graph


def remove_short_fragments(graph, threshold=QUAL.MINIMAL_FRAGMENT_SIZE):
	"""
	Remove fragments which are shorter than 500bp (QUAL.MINIMAL_FRAGMENT_COVERAGE)

	Args:
		graph (decoil.graph.Mutigraph): SV graph

	Returns:
		filtered graph including only fragments with size >= 500 bp
	"""

	# remove low coverage fragments
	fragments = graph.get_fragments()
	toremove = []
	for fid in fragments:
		frag = fragments[fid]
		if frag.len < threshold:
			toremove.append(fid)

	for f in toremove:
		print("Remove fragment due to size", f, fragments[f].coords)
		graph.rewire_fragment_neighbors(f)
		graph.remove_fragment(f)

	return graph


def resolve_self_loops(graph):
	"""
	Resolve fragments having self loops
	"""
	pairs_selfloop = []
	fragments = graph.get_fragments()
	edges = graph.get_edges()
	
	for fid in fragments:
		frag_obj = fragments[fid]
		u = frag_obj.tail
		v = frag_obj.head
		e = frag_obj.edge
		
		for etemp in graph[u][v]:
			# edge connecting head to tail of a fragment other then "SVTYPE = FRAGMENT"
			if etemp != e and edges[etemp][gep.SVTYPE] != gp.FRAGMENT:
				# found self loop
				if edges[etemp][gep.SVTYPE] == gp.DUP:
					print("Found tandem duplication fragid:(tail, head, edge) {}:({},{},{}) {}".format(fid, u, v, e,
					                                                                                   edges[etemp][
						                                                                                   gep.SVTYPE]))
					edges[etemp][gep.REPEAT] = True
				else:
					print("Found other loop type fragid:(tail, head, edge) {}:({},{},{}) {}".format(fid, u, v, e,
					                                                                                edges[etemp][
						                                                                                gep.SVTYPE]))
					edges[etemp][gep.REPEAT] = True
				pairs_selfloop.append((e, u, v))
	
	print("Found loops", pairs_selfloop)
	for e, u, v in pairs_selfloop:
		graph.remove_edge(e, u, v)
	
	return graph


def add_coverage_per_fragment(graph, bigwigfile):
	"""
	Compute the coverage per fragment.

	Args:
		fragment_dict (dict): Fragments with the start and end annotation. Eg. {0:(1, 1999), 1:(2000,3000)}
		graph (decoil.Graph): Structural variants representation
		bigwigfile (pyBigWig.bw): Coverage file

	Returns:
		decail.Graph having for all "fragment" node the weight == coverage
	"""
	bw = pyBigWig.open(bigwigfile)
	
	if bw.isBigWig() == True:
		
		all_chroms = bw.chroms()
		edges = graph.get_edges()
		fragments = graph.get_fragments()
		
		for fid in fragments:
			frag_obj = fragments[fid]
			edge_obj = edges[frag_obj.edge]
			chr, start, stop = frag_obj.chr, frag_obj.start, frag_obj.end
			
			# check for chr nomenclature compatibility
			if chr not in all_chroms:
				raise ValueError(
					"""{} extracted from vcf not present in {}. Please check formatting""".format(str(chr),
					                                                                              bigwigfile))
			
			# for inf values set the length of the chromosome
			if start == math.inf:
				start = all_chroms[chr]
			if stop == math.inf:
				stop = all_chroms[chr]
			# avoid different size of chr
			stop = max(0, min(stop, all_chroms[chr]))
			start = max(0, min(start, all_chroms[chr]))
			
			# compute average across bin
			if start < stop:
				avg_cov = bw.stats(chr, start, stop, type="mean")[0]
			elif start > stop:
				avg_cov = bw.stats(chr, stop, start, type="mean")[0]
			else:
				# equal start and stop
				avg_cov = bw.stats(chr, start - 1, start, type="mean")[0]
			
			if avg_cov > -1:
				# set fragment coverage
				frag_obj.coverage = avg_cov
				# set edge weight as the fragment coverage
				edge_obj.props[gep.WEIGHT] = avg_cov
			# print("Fragment", fid, " coverage: ", fragments[fid].coverage, " edge ", frag_obj.edge, "with coverage",edge_obj.props[gep.WEIGHT])
	else:
		raise ValueError("""{} needs to be a bigWig file""".format(bigwigfile))



def passqc(count_spanning_reads, mean_wgs):
	"""
	Decide if spanning reads are worth to be accounted.

	Arguments:
		count_spanning_reads (int): Number of reads spanning the breakpoint (not split reads)
		mean_wgs (float): WGS mean coverage

	Returns:
		Boolean
	"""
	if count_spanning_reads > 4 and count_spanning_reads > mean_wgs:
		return True
	return False

def add_spatial_edges_new(graph, bamfile, window=500, padd=100):
	"""
	Add edges which connects 2 neighboring fragments in space if supporting reads >= WGS / 2

	Arguments:
		graph (decoil.encode.Multigraph):
		bamfile (str): Path to bam file
		window (int): Left and right window spanning region around breakpoint
		padd (int):
	"""

	log.info("5. Add edges connecting neighboring fragments")

	samfile = pysam.AlignmentFile(bamfile, "rb")
	fragment_intervals = graph.get_fragment_intervals()
	fragments = graph.get_fragments()
	fragments_connected = {}

	frag_keys = list(fragments.keys())

	for i in range(len(frag_keys)-1):

		fid1 = frag_keys[i]
		fid2 = frag_keys[i+1]

		chr1, start1, end1 = fragments[fid1].chr, fragments[fid1].start, fragments[fid1].end
		chr2, start2, end2 = fragments[fid2].chr, fragments[fid2].start, fragments[fid2].end
		breakpoint = -1

		if chr1 == chr2:
			if abs(start1 - end2) < QUAL.DISTANCE:
				breakpoint = start1
			# left = fid2
			# right = fid1
			elif abs(start2 - end1) < QUAL.DISTANCE:
				breakpoint = start2
		# left = fid1
		# right = fid2

		# fragments are neighbors
		if breakpoint > -1:

			count_coverage_through = 0
			for read in samfile.fetch(contig=chr1,
			                          start=max(0, breakpoint - window),
			                          stop=breakpoint + window):

				# filter out secondary alignment and mapped
				if not read.is_mapped or read.is_secondary:
					continue

				# filter out short reads (< 2*window)
				if read.query_length < 2 * window:
					continue

				# keep only reads which go through the breakpoint
				if read.reference_start < (breakpoint - padd) and read.reference_end > (breakpoint + padd):
					count_coverage_through += 1

			if passqc(count_coverage_through,QUAL.MEAN_COVERAGE_WGS):

				# nodes
				left = fragment_intervals.find_head_node(chr1, breakpoint)
				right = fragment_intervals.find_tail_node(chr1, breakpoint)

				settings = {}
				settings["svtype"] = gp.SPATIAL
				settings["weight"] = count_coverage_through
				settings["color"] = "pink"

				conn1 = """{}_{}""".format(left, right)
				conn2 = """{}_{}""".format(right, left)

				# avoid duplicated connection between two fragments
				if left and right and conn1 not in fragments_connected and conn2 not in fragments_connected:
					print("Spatial ->> Chr breakpoint count node_left node_right", chr1, breakpoint,
					      count_coverage_through, left, right)
					graph.add_edge(MultiGraph.generate_id(), left, right, gp.SPATIAL, settings)
					fragments_connected[conn1] = None
					fragments_connected[conn2] = None

	return graph


def run_cleaning(G, outputdir):
	"""
	Cleaning of graph (e.g. remove standalone fragments)
	"""
	os.makedirs(outputdir, exist_ok=True)
	log.info("4.1 Remove lowcoverage fragments")
	G = remove_lowcoverage_fragments(G, threshold=QUAL.MINIMAL_FRAGMENT_COVERAGE)

	log.info("4.1 Remove lowcoverage fragments")
	G = remove_short_fragments(G, threshold=QUAL.MINIMAL_FRAGMENT_SIZE)

	log.info("4.3 Remove standalone fragments")
	G = remove_standalone_fragments(G)
	
	# log.info("5. Balance graph")
	# G = encode.balance(G)
	
	# log.info("5. Remove self loops")
	# G = encode.resolve_self_loops(G)
	
	# log.info("5.1 Debug fragments")
	
	return G
