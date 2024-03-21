"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 5:45 PM 9/19/22

Search cycles in a graph business logic.
"""
import math
import numpy as np
from copy import deepcopy
from collections import defaultdict
from decoil.utils import PATH
from decoil.validate import compare as compare


def is_overlap(p1, p2):
	"""
	Check if the end of fragment p1 overlaps with the start fragment p2
	"""
	if p1[PATH.CHR] == p2[PATH.CHR] and p1[PATH.STRAND] == "+" and p1[PATH.STRAND] == p2[PATH.STRAND] and (
		abs(p1[PATH.END] - p2[PATH.START]) < 10):
		return True

	if p1[PATH.CHR] == p2[PATH.CHR] and p1[PATH.STRAND] == "-" and p1[PATH.STRAND] == p2[PATH.STRAND] and (
		abs(p1[PATH.START] - p2[PATH.END]) < 10):
		return True

	return False

def do_collapse(p1, p2):
	"""
	Merge p2 into p1
	"""
	frag = """{}.{}""".format(str(p1[PATH.FRAG]), str(p2[PATH.FRAG]))
	cov = (p1[PATH.COV] + p2[PATH.COV]) / 2

	if p1[PATH.STRAND] == "+":
		start = p1[PATH.START]
		end = p2[PATH.END]
	else:
		start = p2[PATH.START]
		end = p1[PATH.END]

	return (p1[PATH.CHR], start, end, p1[PATH.STRAND], frag, p1[PATH.CIRCID], cov)


def rotate_toleft(path):
	if len(path) > 1:
		newpath = path[1:]
		newpath.append(path[0])
		return newpath
	return path


def rotate_toright(path):
	if len(path) > 1:
		newpath = [path[-1]] + path[:len(path) - 1]
		return newpath
	return path


def rotate_tominimal(path):
	"""
	Rotate path until the left most fragment is also genomically left most.
	Assumes frag id are integers and relies on the property that min(frag)==genomically left most fragment

	Arguments:
		paths (list): List fragments (cycle)
	"""
	# important to consider first hit as being the smaller one
	#                       ---------
	#
	lpath  = [abs(p) for p in path]
	min = math.inf
	idx = -1
	for i, p in enumerate(lpath):
		if min > p:
			min = p
			idx = i
	return path[idx:] + path[:idx]


def collapse_fragments(path):
	"""
	Collapse together consecutive fragments for a cycle path.

	[('chr2', 15000001, 15102001, '-', 5, 'circ_0_1', 45.34656862745098),
	('chr2', 15310000, 15390001, '-', 17, 'circ_0_1', 82.17237284533944),
	('chr2', 15102001, 15120002, '-', 9, 'circ_0_1', 86.67751791567134)]

	transform to

	[('chr2', 15310000, 15390001, '-', 17, 'circ_0_1', 82.17237284533944),
	('chr2', 15000001, 15120002, '-', '9.5', 'circ_0_1', 66.01204327156117)]

	"""
	if len(path) == 1:
		return path
	else:
		newp = []
		current_fragment = path[0]
		for i in range(1, len(path)):
			# same chromosome and same strand and end(current_fragment) overlaps with start(p[i+1])
			if is_overlap(current_fragment, path[i]):
				# merge current_fragment with p[i]
				current_fragment = do_collapse(current_fragment, path[i])
			else:
				newp.append(current_fragment)
				current_fragment = path[i]

		# if there is an overlap between first and last fragment
		if is_overlap(current_fragment, path[0]):
			current_fragment = do_collapse(current_fragment, path[0])
			# skip first element in the list because it was fused to the last
			newp.append(current_fragment)
			newp = newp[1:]
		else:
			newp.append(current_fragment)

		return newp

def flip_strands(paths):
	newpaths = []
	for p in paths:
		count_plus = len([e[PATH.STRAND] for e in p if e[PATH.STRAND] == "+"])
		count_minus = len([e[PATH.STRAND] for e in p if e[PATH.STRAND] == "-"])

		# more positive strands
		if count_plus > count_minus:
			newpaths.append(p)
		else:
			newp = []
			for e in p[::-1]:
				strand = "+"
				if e[PATH.STRAND] == "+":
					strand = "-"
				newp.append((e[PATH.CHR],
				             e[PATH.START],
				             e[PATH.END],
				             strand,
				             e[PATH.FRAG],
				             e[PATH.CIRCID],
				             e[PATH.COV]
				             ))
			newpaths.append(newp)

	return newpaths

def flip_strand(path):
	"""
	Flip strands of cycles (in case more inverted fragments then forward fragments)
	"""
	count_plus = len([f for f in path if f >= 0])

	# more inverted fragments
	if count_plus < (len(path) - count_plus):
		return [-f for f in path[::-1]]
	else:
		return path

def format_paths(paths):
	"""
	Format circular path (flip strands and rotate path to have left most fragment
	"""
	# flip strands (more positive then reverse)
	n1 = flip_strands(paths)

	# rotate to minimal fragment (convention for circle start)
	n2 = rotate_tominimal(n1)

	return n2


def simplify_paths(paths):
	"""
	Simplify paths (correct overfragmentation of the paths, flip strands)
	"""

	# flip strands (more positive then reverse)
	n1 = flip_strands(paths)

	# rotate to minimal fragment (convention for circle start)
	n2 = rotate_tominimal(n1)

	# collapse fragments which are consecutive in path and in genome
	n3 = collapse_fragments(n2)

	return n3

def transform_fid2genomicpos(candidates, graph):
	"""
	Expand information per cycle from [5,17, 9] to

	[('chr2', 15000001, 15102001, '-', 5, 'circ_0_1', 45.34656862745098),
	('chr2', 15310000, 15390001, '-', 17, 'circ_0_1', 82.17237284533944),
	('chr2', 15102001, 15120002, '-', 9, 'circ_0_1', 86.67751791567134)]
	"""
	fragments = graph.get_fragments()
	expanded_candidates = defaultdict(dict)

	for c in candidates:
		path = []
		for fid in candidates[c]["conf"]:
			fobj = fragments[abs(fid)]
			strand = "+" if fid >= 0 else "-"
			chr, start, stop, coverage = fobj.chr, fobj.start, fobj.end, fobj.coverage
			path.append((chr,
			             start,
			             stop,
			             strand,
			             abs(fid),
			             c,
			             coverage))
		expanded_candidates[c]["score"] = candidates[c]["score"]
		expanded_candidates[c]["conf"] = path

	return expanded_candidates

def simplify_likely_candidates(candidates):
	"""
	Join together fragments if next to each other in the genome.
	These candidates come from the graph.metahash which are already canonical conformation
	"""
	for circ in candidates:
		# # flip strands (more positive then reverse)
		# n1 = flip_strands([candidates[circ]["conf"]])
		#
		# # rotate to minimal fragment (convention for circle start)
		# n2 = rotate_tominimal(n1)

		# collapse fragments which are consecutive in path and in genome
		n3 = collapse_fragments(candidates[circ]["conf"])

		# reassign paths
		candidates[circ]["conf"] = n3

	return candidates



def format_combined_paths(paths):
	"""
	Format paths (correct overfragmentation of the paths, flip strands)

	Arguments:
		paths (dict):
	"""
	circ_id = [circ for circ in paths]
	for circ in circ_id:
		# flip strands (more positive then reverse)
		n1 = flip_strands([paths[circ]])

		# rotate to minimal fragment (convention for circle start)
		n2 = rotate_tominimal(n1)

		# reassign paths
		paths[circ] = n2[0]

	return paths


def deduplicate_combined_paths(paths):
	"""
	Remove duplicates from paths
	Consider equal paths only if lists are identical (DO NOT take into account circularity
	"""

	deduplicate = {}

	for key in paths:
		if len(deduplicate) == 0:
			deduplicate[key] = paths[key]
		else:
			p1 = paths[key]
			duplicate = False
			for key2 in deduplicate:
				p2 = deduplicate[key2]
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
				deduplicate[key] = p1

	return deduplicate


def insert_edges_sorted(list_edges, e, data):
	"""
	Insert edge following criteria
		- edges weight (descendent)
		- for same weight consider sorting:
			- connecting fragment coverage (descending)
			- if 2 edges connect 2 same nodes, with same weight then priority DEL, DUP, INV, INS, BND, SPATIAL

	Arguments:
		list_edges (list): List of tuples
		e (int): Edge id
		data (dict): Information about the edge (see __collect_edge_info__)
					data = {"rank1": 150,
							"rank2": 50,
							"rank3": 3,
							"v": 3}
	"""
	# insertion index
	idx = -1

	# no elements in the list
	if len(list_edges) == 0:
		return list_edges.append((e, data))

	# find insertion point
	for i, e_temp in enumerate(list_edges):
		# find insertion spot using max coverage
		if e_temp[1]["rank1"] <= data["rank1"]:
			idx = i
			break

	# check if last position
	if idx == -1:
		list_edges.append((e, data))
		return list_edges

	# insertion point found due to different rank1
	if list_edges[i][1]["rank1"] < data["rank1"]:
		list_edges.insert(i, (e, data))
		return list_edges

	# for equal rank1 find insertion point based on rank2
	j = i
	while j < len(list_edges) and \
		list_edges[j][1]["rank1"] == data["rank1"] and \
		list_edges[j][1]["rank2"] > data["rank2"]:
		j += 1

	# equal rank1, differnt rank2 but reaching last element
	if j == len(list_edges):
		list_edges.append((e, data))
		return list_edges

	# equal rank1 and different rank2
	a1 = list_edges[j][1]["rank2"]
	a2 = data["rank2"]
	if a1 < a2:
		list_edges.insert(j, (e, data))
		return list_edges

	# equal rank1 and different rank2
	if a1 > a2:
		if j + 1 < len(list_edges):
			list_edges.insert(j + 1, (e, data))
		else:
			list_edges.append((e, data))
		return list_edges

	# equal rank1, equal rank2, different rank3
	if data["rank3"] < list_edges[j][1]["rank3"]:  # the smaller to left
		list_edges.insert(j, (e, data))
		return list_edges

	# all equal
	if (j + 1) < len(list_edges):
		list_edges.insert(j + 1, (e, data))
	else:
		list_edges.append((e, data))
	return list_edges


def get_hash(p):
	"""
	Compute hash for a list.

	Arguments:
		p (list): List of fragments ids

	Returns:
		hash value
	"""
	return hash(tuple(p))


def get_metahash(chain_fragments):
	"""
	Compute metahash of a circular structure (consider all rotations and reverse complement).

	Arguments:
		chain_fragments (list): Chain of fragments ids describing one cycle

	Returns:
		hash value
	"""
	forward = chain_fragments
	reverse = [-f for f in chain_fragments[::-1]]

	# initialize hash
	metahash = hash(tuple(forward)) ^ hash(tuple(reverse))

	for conf in [forward, reverse]:
		for i in range(1, len(conf)):
			rotation = conf[i:] + conf[:i]
			metahash ^= hash(tuple(rotation))
	return metahash

def convert2canonical_path(path):
	"""
	Rotate path until the left most fragment is also genomically left most.
	Assumes frag id are integers and relies on the property that min(frag)==genomically left most fragment

	Arguments:
		paths (list): List of fragments representing one cycle

	"""

	# flip strands (more positive then reverse)
	n1 = flip_strand(path)

	# rotate to minimal fragment (convention for circle start)
	n2 = rotate_tominimal(n1)

	return n2

def filter_reconstructions(dict_cycle_candidates, graph, threshold = 0.7):
	"""
	For similar structures keep the one with highest score
	"""
	id_circles = list(dict_cycle_candidates.keys())
	id_circles_visited = np.zeros(len(id_circles))
	dict_fragments = graph.get_fragments()

	for i in range(0, len(id_circles)):
		for j in range(i, len(id_circles)):
			if i != j:
				mcid = id_circles[i]
				cid = id_circles[j]
				jscore, s1, s2 = compare.compute_jaccard(dict_cycle_candidates[mcid]["conf"],
														 dict_cycle_candidates[cid]["conf"], graph)
				if jscore > threshold:
					if dict_cycle_candidates[mcid]["score"] < dict_cycle_candidates[cid]["score"]:
						# set circle as similar / to remove
						id_circles_visited[i] = 1
					else:
						id_circles_visited[j] = 1

	to_remove = []
	for i in range(0,len(id_circles_visited)):
		if id_circles_visited[i] == 1:
			to_remove.append(id_circles[i])

	for k in to_remove:
		del dict_cycle_candidates[k]

	return dict_cycle_candidates
