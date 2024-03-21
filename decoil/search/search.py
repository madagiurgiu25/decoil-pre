"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 5:30 PM 9/19/22
"""

import logging
import os
import sys
from collections import defaultdict
from copy import deepcopy
from itertools import combinations

import decoil.utils as utils
from decoil.search import cycles
from decoil.utils import GRAPH_NODE_PROP as gnp
from decoil.utils import GRAPH_PROP as gp
from decoil.validate import compare as compare

log = logging.getLogger('reconstruct.search')
log.propagate = False
log.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
log.addHandler(handler)


def add_fragment(v1, v2):
	chr1, pos1, part1 = v1.split(utils.SEPARATOR)
	chr2, pos2, part2 = v2.split(utils.SEPARATOR)
	
	# the chr should be the same
	if chr1 != chr2:
		return None, None, None, None
	chr = chr1
	
	if part1 == 't' and part2 == 't':
		return None, None, None, None
	
	if part1 == 'h' and part2 == 'h':
		return None, None, None, None
	
	# test forward
	if part1 == 't' and part2 == 'h':
		strand = "+"
	else:
		strand = "-"
	start = min(int(pos1), int(pos2))
	end = max(int(pos1), int(pos2))
	
	return chr, start, end, strand


def chain_fragments(paths, graph):
	"""
	Chain fragments including genomic information for the paths
	
	Returns:
		List of paths
		[ [(chr, start, stop, strand, frag, circid),
		   (chr, start, stop, strand, frag, circid)]
		]
	"""
	paths_long = []
	for j, p in enumerate(paths):
		new_path_long = []
		# read node pairs
		for i in range(len(p) - 1):
			(n1, e1) = p[i]
			(n2, e2) = p[i + 1]
			
			if graph.edges[e2].svtype == gp.FRAGMENT:
				
				_chr = graph.nodes[n1].props[gnp.CHR]
				_pos1 = graph.nodes[n1].props[gnp.POS]
				_pos2 = graph.nodes[n2].props[gnp.POS]
				_frag = graph.nodes[n2].parent_fragment
				_cov = graph.edges[e2].weight
				
				if _pos1 < _pos2:
					strand = "+"
				else:
					strand = "-"
					tmp = _pos1
					_pos1 = _pos2
					_pos2 = tmp
				
				new_path_long.append((_chr,
				                      _pos1,
				                      _pos2,
				                      strand,
				                      _frag,
				                      "circ_{}".format(j),
				                      _cov))
		paths_long.append(new_path_long)
	
	return paths_long


def format_and_check(path, graph, circid="default"):
	"""
	Format the path as a list of fragments (with orientation)
	"""
	curated_path = []
	
	start_fragment = False
	start_node = path[0][0]
	end_node = path[-1][1]
	state = -1  # 0 - link, 1 - fragment
	
	for i, (v1, v2, key) in enumerate(path):
		# initial
		if i == 0:
			if graph[v1][v2][key]["svtype"] == 'fragment':
				# mark that first edge is a fragment
				start_fragment = True
				state = 1
			else:
				state = 0
		else:
			# coming from a link to a fragment
			if state == 0 and graph[v1][v2][key]["svtype"] == 'fragment':
				state = 1
				chr, start, stop, strand = add_fragment(v1, v2)
				if chr is not None:
					curated_path.append((chr, start, stop, strand, graph[v1][v2][key]["frag"], circid))
			else:
				# coming from fragment to link
				if state == 1 and graph[v1][v2][key]["svtype"] != 'fragment':
					state = 0
				else:
					log.error("The path is not clean")
					exit(1)
	
	# check circularity
	if start_fragment == False and state == 1 and start_node == end_node:
		return curated_path
	if start_fragment == True and state == 0 and start_node == end_node:
		return curated_path
	
	# if something went wrong
	return []


def write_rawpath(path, file):
	with open(file, "w") as out:
		for item in path:
			out.write("{}\t{}\n".format(item[0], item[1]))


def get_longest_path(found_paths, gcomp):
	size = 0
	index = -1
	for i, p in enumerate(found_paths):
		count_frag = 0
		for (u, v, k) in p:
			if gcomp[u][v][k]["svtype"] == "fragment":
				count_frag += 1
		if size < count_frag:
			index = i
	if index != -1:
		return deepcopy(found_paths[index])
	else:
		return []


def deduplicate(parent_list, new_entries):
	"""
	Attach to parent_list only circular paths not existant in new_entries
	Consider equal paths only if lists are identical (DO NOT take into account circularity
	"""
	tokeep = []

	if len(parent_list) == 0:
		return new_entries

	for p1 in new_entries:
		duplicate = False
		for p2 in parent_list:
			if len(p1) != len(p2):
				continue
			match = True
			for i in range(len(p1)):
				# check if path contains same nodes and edges
				if i == 0 and p1[i][0] != p2[i][0]:
					match = False
					break
				elif i != 0 and (p1[i][0] != p2[i][0] or p1[i][1] != p2[i][1]):
					match = False
					break

			# path already exists
			if match == True:
				duplicate = True
				break

		# no duplication found
		if not duplicate:
			tokeep.append(p1)

	return parent_list + tokeep

def get_start_nodes(allnodes):
	"""
	Get all nodes from which you can start searching
	"""
	nodes = allnodes
	narr = []
	for n in nodes:
		if nodes[n].part == gp.TAIL:
			narr.append(n)
	return narr

def is_fragment_consumed(graph, f, fobj, explored_paths):
	"""
	Check if fragment already fully explored
	"""
	tail, head  = fobj.tail, fobj.head

	# loose end not worth exploring for circular path
	if graph.degree(tail) == 0 or graph.degree(head) == 0:
		return True

	# in = out == 1 means
	if graph.degree(tail) == 1 and graph.degree(head) == 1:
		for p in explored_paths:
			if f in p or -f in p:
				return True
	return False

def find_path(graph):
	"""
	Find circular path by DFS tree
	"""
	toreturn = []
	found_paths = []
	allfragments  = graph.get_fragments()
	for f in allfragments:
		if not is_fragment_consumed(graph,f, allfragments[f],found_paths):
			found_paths = graph.find_simple_circles3(allfragments[f].tail)
		# found_paths = graph.find_simple_circles(v)
		# found_paths = graph.find_simple_circles2(v)
			if len(found_paths) > 0:
				toreturn += found_paths
	return toreturn


def run_searching(G, outputdir, name):
	"""
	Search cycles in the graph by a DFS starting from every node
	"""
	os.makedirs(outputdir, exist_ok=True)
	# pathfileraw = os.path.join(outputdir, "reconstruct_path.txt")
	
	log.info("6. Find path")
	findings = find_path(G)  # compute cycle search across all graph connected components
	return findings


def run_searching_v2(dfs):
	"""
	Search for all circular paths.

	Returns:
		list of circular paths (fragments chains, strand orientation aware)
	"""
	log.info("5. Compute simple circles")
	# compute all cycles
	dfs.reset_cycles()
	dfs.compute_all_cycles(dfs.root, None, None, None, None, [])
	print(dfs.cycles)
	
	# chain fragments
	dfs.reset_fragments()
	dfs.cycles_to_fragments()
	
	return dfs.fragments


def lookup_unique(p , graph):
	"""
	Check if the cycle already indexed (consider rotations and reverse complement).
	"""
	h = cycles.get_metahash(p)
	if h in graph.get_metahash_circles():
		return False
	return True

def lookup_unique_canonical(p , graph):
	"""
	Check if the cycle already indexed.

	Arguments:
		p (list): Cycle given in the canonical form
	"""
	h = cycles.get_hash(p)
	if h in graph.get_hash_circles():
		return False
	return True

def lookup_different_enough(c1, list_cycles, graph, threshold = 0.7):

	for c2 in list_cycles:
		jscore, _, _ = compare.compute_jaccard_unique_fragments(c1,c2, graph)
		if jscore > threshold:
			return False
	return True


def run_combine_paths(paths, graph):
	"""
	Combine simple circle paths.
	"""
	log.info("6. Combine simple circles")
	paths_starts = defaultdict(list)
	derived_circles = []

	# Group cycles starting with same fragment
	for path in paths:
		paths_starts[path[0]].append(path)

	# Filter out highly similar cycles
	for p in paths_starts:
		paths_starts[p] = compare.filter_out_similar_cycles(paths_starts[p], graph)


	# Compute all permutations/combinations for paths starting with same id
	for p in paths_starts:

		if len(paths_starts[p]) >= 2:
			# find all combinations of paths starting with same fragment
			# total circles to combine
			total_len = len(paths_starts[p])
			index_list = range(total_len)

			# find comb_n,1,comb_n,2
			for i in range(2, total_len + 1):

				# compute all combinations i from n
				comb = list(combinations(index_list, i))
				for ctemp in comb:
					newpath = []
					# extract the positions, e.g. 0,1,2
					for j in ctemp:
						_circle = paths_starts[p][j]
						newpath = newpath + _circle

					newpath_canonical = cycles.convert2canonical_path(newpath)
					if lookup_unique_canonical(newpath_canonical, graph) and lookup_different_enough(newpath, derived_circles, graph):
						derived_circles.append(newpath)
						graph.add_to_hash(newpath)
