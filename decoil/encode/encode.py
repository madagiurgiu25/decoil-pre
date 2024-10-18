"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 9:30 AM 9/19/22
"""
import os
import numpy
import logging
import math
import sys
from collections import defaultdict

import vcfpy
from vcfpy import InvalidRecordException

import decoil.encode.operations as operations
import decoil.utils as utils
from decoil.encode.graph import MultiGraph
from decoil.utils import FRAGMENT as fg
from decoil.utils import GRAPH_EDGE_PROP as gep
from decoil.utils import GRAPH_NODE_PROP as gnp
from decoil.utils import GRAPH_PROP as gp
from decoil.utils import VCF_PROP as vp

log = logging.getLogger('reconstruct.encode')
log.propagate = False
log.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)

def cleanvcf(vcfin, vcfout, multi=False):
	"""
	Remove corrupted entries in VCF file (generated by Sniffles)
	"""
	with open(vcfout, "w") as g:
		with open(vcfin, "r") as f:
			for l in f:
				if l.startswith("#"):
					g.write(l)
				else:
					if multi == False and len(l.strip().split("\t")) == 10:
						g.write(l)
					elif multi == False and len(l.strip().split("\t")) != 10:
						log.warning("Corrupted VCF line " + l)
					else:
						g.write(l)

def transform_record(record, svcaller):
	"""
	Convert required fields from Sniffles2 or others to Sniffles
	"""
	if svcaller == vp.CUTESV:
		# add STRANDS instead of STRAND or add an artificial strand encoding for INS
		record.INFO[vp.STRANDS] = [record.INFO[vp.STRAND] if vp.STRAND in record.INFO else "++"]

		# add CHR2 and POS2
		record.INFO[vp.CHR2] = record.CHROM if record.INFO[vp.SVTYPE] not in [vp.BND, vp.TRA] else record.ALT[0].mate_chrom
		record.INFO[vp.END] = record.INFO[vp.END] if record.INFO[vp.SVTYPE] not in [vp.BND, vp.TRA] else record.ALT[0].mate_pos

	if svcaller == vp.SNIFFLES2:
		# add STRANDS instead of STRAND
		record.INFO[vp.STRANDS] = [record.INFO[vp.STRAND]]

		# add CHR2 and POS2
		record.INFO[vp.CHR2] = record.CHROM if record.INFO[vp.SVTYPE] not in [vp.BND, vp.TRA] else record.ALT[0].mate_chrom
		record.INFO[vp.END] = record.INFO[vp.END] if record.INFO[vp.SVTYPE] not in [vp.BND, vp.TRA] else record.ALT[0].mate_pos

		if record.calls[0].data[vp.GT] == "./.":
			record.calls[0].data[vp.DV] = int(record.INFO[vp.SUPPORT])
			record.calls[0].data[vp.DR] = int(record.INFO[vp.COVERAGE][0]) if record.INFO[vp.COVERAGE][0] != 'None' else 0
			record.calls[0].data[vp.AF] = (record.calls[0].data[vp.DV]/(record.calls[0].data[vp.DV]+record.calls[0].data[vp.DR]+1))

			if record.calls[0].data[vp.AF] <= 0.3:
				record.calls[0].data[vp.GT] = "0/0"
			elif record.calls[0].data[vp.AF] <= 0.7:
				record.calls[0].data[vp.GT] = "0/1"
			else:
				record.calls[0].data[vp.GT] = "1/1"

	return record

def readvcf(vcffile, outputdir, multi=False, svcaller=vp.SNIFFLES):
	"""
	Read vcf file and store all breakpoints.

	Returns:
		(breakpoints, sv info)
	"""
	# breakpoints info
	collection_breakpoints = defaultdict(list)  # chr2: [pos1, pos2, pos3]
	# sv info
	svinfo = defaultdict(list)  # "chr2:pos": [(key1, SVTYPE, coverage), (key2, SVTYPE2, coverage)]

	# check for corrupted entries
	vcffile_clean = os.path.join(outputdir,"clean.vcf")
	cleanvcf(vcffile,vcffile_clean, multi=multi)

	# parse vcf
	reader = vcfpy.Reader.from_path(vcffile_clean)
	count=0
	try:
		for record in reader:
			print(record)
			if svcaller in [vp.SNIFFLES2, vp.CUTESV]:
				record = transform_record(record, svcaller)

			id = record.ID[0]  # primary id
			count+=1
			svtype = record.INFO[vp.SVTYPE]
			strand = record.INFO[vp.STRANDS][0]

			# single sample vcf
			if multi == False:
				dv = record.calls[0].data.get(vp.DV)
				dr = record.calls[0].data.get(vp.DR)
				gt = record.calls[0].data.get(vp.GT)
			# multi-sample vcf
			else:
				dv_arr = [record.calls[i].data['DR'][1] for i in range(len(record.calls))]
				dr_arr = [record.calls[i].data['DR'][0] for i in range(len(record.calls))]
				gt_arr = [record.calls[i].data['GT'] for i in range(len(record.calls))]
				max_index = numpy.argmax(numpy.array(dv_arr))
				dv = dv_arr[max_index]
				dr = dr_arr[max_index]
				gt = gt_arr[max_index]

			# check quality filter
			if not operations.pass_filter(record, multi=multi):
				# skip record as this is good
				continue

			if svtype in [vp.DUP, vp.DEL, vp.INV, vp.BND, vp.INVDUP, vp.DUPINV, vp.TRA]:
				chr1 = str(record.CHROM)
				pos1 = str(record.POS)
				
				chr2 = str(record.INFO[vp.CHR2])
				pos2 = str(record.INFO[vp.END])

				print(chr1, pos1, chr2, pos2)
			
			# elif svtype in [vp.INS, vp.DUPINS]:  # ! sniffles encodes SVTYPE==INS for both INS and INDELS!
			# 	chr1 = str(record.CHROM)
			# 	pos1 = str(record.POS)
			#
			# 	chr2 = str(record.INFO[vp.CHR2])
			# 	# len of insertion
			# 	inslen = len(record.INFO[vp.SEQ][0]) if vp.SEQ in record.INFO else 0
			# 	inslen += len(record.ALT[0].value) if isinstance(record.ALT[0], vcfpy.Substitution) else 0
			# 	pos2 = str(record.POS + inslen)
			
			# DV - coverage alternative allele, DR - coverage reference allele
			svinfo['{}{}{}'.format(chr1, utils.SEPARATOR, pos1)].append((id, chr2, pos2, svtype, dv, dr, gt, strand))
			collection_breakpoints[chr1].append(int(pos1))
			collection_breakpoints[chr2].append(int(pos2))
	
	except InvalidRecordException:
		log.warning("VCF file was cropped! Most probably because Sniffle VCF is corrupted. This will generate incomplete/incorrect reconstruction.")
		pass

	print("Total number of entries in vcf", count)
	print("SV kept for graph generation", len(svinfo))
	# sort breakpoints
	for key in collection_breakpoints:
		collection_breakpoints[key] = sorted(list(set(collection_breakpoints[key])))
	
	# log.info("1. 1. clean breakpoints")
	# log.info(collection_breakpoints)
	# refine breakpoints
	collection_breakpoints, svinfo = operations.merge_near_breakpoints(collection_breakpoints, svinfo)
	log.info("1. 2. After cleaning breakpoints")
	# print(collection_breakpoints)
	# print(svinfo)
	
	return collection_breakpoints, svinfo


def insert_fragment(graph, chr, pos, last_pos, fragment_count):
	"""
	Add "fragment" node in the graph.
	"""
	
	# chr1:12345:t - from tail
	# chr1:12345:h - from head
	if pos == math.inf:
		pos = utils.POS.MAX
	
	id1 = MultiGraph.generate_id()
	id2 = MultiGraph.generate_id()
	id3 = MultiGraph.generate_id()
	
	props1 = {gnp.CHR: chr, gnp.POS: last_pos, gnp.COLOR: "black"}
	props2 = {gnp.CHR: chr, gnp.POS: pos, gnp.COLOR: "black"}
	
	graph.add_node(id1, gp.TAIL, fragment_count, props1)
	graph.add_node(id2, gp.HEAD, fragment_count, props2)
	
	frag_len = (int(pos) - int(last_pos))
	props3 = {gep.CHR: chr,
	          gep.START: last_pos,
	          gep.END: pos,
	          gep.COLOR: "gray",
	          gep.SVTYPE: gp.FRAGMENT,
	          gep.FRAGMENT_LEN: frag_len,
	          gep.FRAGMENT_NAME: fragment_count,
	          gep.WEIGHT: math.inf}
	graph.add_edge(id3, id1, id2, gp.FRAGMENT, props3)
	
	graph.add_fragment(fragment_count, id1, id2, id3, {fg.CHR: chr,
	                                                   fg.START: last_pos,
	                                                   fg.END: pos,
	                                                   fg.COVERAGE: math.inf
	                                                   })
	graph.get_fragment_intervals().add_fragment(chr, last_pos, pos, {fg.FRAGID: fragment_count,
	                                                                 fg.EDGEID: id3,
	                                                                 fg.NODE_TAIL: id1,
	                                                                 fg.NODE_HEAD: id2
	                                                                 })
	return graph


# def connect_two_fragments(graph, left_fragment_id, right_fragment_id, edge_id):
#
# 	if left_fragment_id != -1 and right_fragment_id != -1:
# 		left_node_id = graph.get_fragments()[left_fragment_id][gp.HEAD]
# 		right_node_id = graph.get_fragments()[right_fragment_id][gp.TAIL]
#
# 		props = {}
# 		props["svtype"] = gp.SPATIAL
# 		props["weight"] = -1  # initial value
# 		props["color"] = "pink"
#
# 		print("connect fragments ",left_node_id, right_node_id)
# 		graph.add_edge(edge_id, left_node_id, right_node_id, gp.SPATIAL, props)


def fragment_genome(breakpoints, bigwigfile, bamfile):
	"""
	Fragment genome using all the SV breakpoints and encode these into a graph
	Nodes are breakpoints, Edges {SV, Fragment}
	Every "Fragment" Edge has a "head" and "tail" node to mark the orientation of the fragment

	Returns:
		decoil.MultiGraph
	"""
	G = MultiGraph()
	
	# for every chromosome
	for j, chr in enumerate(breakpoints):
		last_pos = 0
		
		# for every SV breakpoint
		for i, pos in enumerate(breakpoints[chr]):
			# insert fragment as node in graph
			current_id = MultiGraph.generate_id()
			G = insert_fragment(G, chr, pos, last_pos, current_id)
			last_pos = pos
		
		# last fragment as node
		G = insert_fragment(G, chr, math.inf, last_pos, MultiGraph.generate_id())
	
	# update all the node pairs
	G.update_pairs()
	
	# add coverage
	operations.add_coverage_per_fragment(G, bigwigfile)
	
	return G


def is_hom(gt, dr, dv):
	"""
	Check the genotype of the variants.
	1/1 - hom
	0/1 - het
	"""
	if dv / (dr + dv) > 0.7:
		return True
	return False


def add_neighbors(graph, node, settings, breakpoints, svtype):
	"""
	Connect spatially close fragments.
	Use frag variable to find the neighbors
	"""
	
	# check if head or tail
	nodes = graph.get_nodes()
	state = nodes[node].part
	
	chr = nodes[node].props[gnp.CHR]
	pos = nodes[node].props[gnp.POS]
	
	fragment_intervals = graph.get_fragment_intervals()
	
	settings["svtype"] = gp.SPATIAL
	settings["weight"] = settings["dr"]
	settings["color"] = "pink"
	
	if state == gp.HEAD:
		# connect to tail of the right neighbor
		nextnode = fragment_intervals.find_tail_node(chr, pos)
	# nextnode = """{}{}{}{}t""".format(chr, utils.SEPARATOR, pos, utils.SEPARATOR)
	else:
		# connect to head of left neighbor
		nextnode = fragment_intervals.find_head_node(chr, pos)
	# nextnode = """{}{}{}{}h""".format(chr, utils.SEPARATOR, pos, utils.SEPARATOR)
	
	# print("add neighbor to node ", node, " the ", nextnode)
	if nextnode:
		graph.add_edge(MultiGraph.generate_id(), node, nextnode, gp.SPATIAL, settings)
	
	return graph


def addsv(graph, svinfo):
	"""
	Encode the SV's as edges in the graph.
	"""
	
	fragments_intervals = graph.get_fragment_intervals()
	
	for key in svinfo:
		for (id, chr2, pos2, svtype, dv, dr, gt, strand) in svinfo[key]:
			chr1, pos1 = key.split(utils.SEPARATOR)
			
			data = {"dv": dv,
			        "dr": dr,
			        "total_cov": dv + dr,
			        "len": 10,
			        "color": vp.SV_COLORS[svtype],
			        "weight": dv}
			
			svtype_code = -1
			if svtype in [vp.BND, vp.TRA]:
				
				if strand == "-+":
					# -+
					# tail head
					v1 = fragments_intervals.find_tail_node(chr1, pos1)
					v2 = fragments_intervals.find_head_node(chr2, pos2)
					svtype_code = gp.BND

				elif strand == "+-":
					# '+-'
					# head tail
					v1 = fragments_intervals.find_head_node(chr1, pos1)
					v2 = fragments_intervals.find_tail_node(chr2, pos2)
					svtype_code = gp.BND

				elif strand == 	"++":
					# '++'
					# head head
					v1 = fragments_intervals.find_head_node(chr1, pos1)
					v2 = fragments_intervals.find_head_node(chr2, pos2)
					svtype_code = gp.BND
				else:
					# '--'
					# tail tail
					v1 = fragments_intervals.find_tail_node(chr1, pos1)
					v2 = fragments_intervals.find_tail_node(chr2, pos2)
					svtype_code = gp.BND

				
				graph.add_edge(MultiGraph.generate_id(), v1, v2, svtype_code, data)
			
			elif svtype in [vp.INS, vp.DEL]:
				# +-
				# head tail
				v1 = fragments_intervals.find_head_node(chr1, pos1)
				v2 = fragments_intervals.find_tail_node(chr2, pos2)
				# print("INS/DEL +-", v1, v2)
				if svtype == vp.INS:
					svtype_code = gp.INS
				else:
					svtype_code = gp.DEL
				
				graph.add_edge(MultiGraph.generate_id(), v1, v2, svtype_code, data)
			
			elif svtype == vp.INV:
				# ++
				# head head
				v1 = fragments_intervals.find_head_node(chr1, pos1)
				v2 = fragments_intervals.find_head_node(chr2, pos2)
				svtype_code = gp.INV

				graph.add_edge(MultiGraph.generate_id(), v1, v2, svtype_code, data)
				
				# --
				# tail tail
				v1 = fragments_intervals.find_tail_node(chr1, pos1)
				v2 = fragments_intervals.find_tail_node(chr2, pos2)
				svtype_code = gp.INV

				graph.add_edge(MultiGraph.generate_id(), v1, v2, svtype_code, data)
			
			elif svtype == vp.DUP:
				# -+
				# tail head
				v1 = fragments_intervals.find_tail_node(chr1, pos1)
				v2 = fragments_intervals.find_head_node(chr2, pos2)
				svtype_code = gp.DUP
				
				graph.add_edge(MultiGraph.generate_id(), v1, v2, svtype_code, data)


		# if SV heterozygous insert spatial relations
		# if is_hom(gt, dr, dv) == False:
		# 	print("SV id:", id, svtype, "is het!!, add neighbors for ", v1, v2)
		# graph = add_neighbors(graph, v1, data, breakpoints, svtype)
		# graph = add_neighbors(graph, v2, data, breakpoints, svtype)
	
	return graph


def run_encoding(vcffile, bigwigfile, bamfile, outputdir,  multi=False, svcaller="sniffles"):
	"""
	Encode the SV information in graph representation
	"""
	log.info("0. Start processing ")
	log.info(vcffile)
	log.info("1. Get all breakpoints and parse SV information")
	breakpoints, svinfo = readvcf(vcffile, outputdir, svcaller=svcaller, multi=multi)
	
	log.info("2. Fragment the genome")
	G = fragment_genome(breakpoints, bigwigfile, bamfile)
	
	log.info("3. Encode SV;s in graph")
	G = addsv(G, svinfo)

	return G
