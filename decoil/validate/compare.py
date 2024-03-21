"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 09:30 AM 02/03/23

Definition and methods of the underlying graph.
"""
import os
import pprint
import numpy as np
import pandas as pd
from collections import defaultdict

import subprocess
from io import StringIO
from intervaltree import Interval, IntervalTree
import pybedtools
import uuid

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from decoil.output import metrics
from decoil.utils import BED_PROP as bp

class Prop():

	def __init__(self, fid, chr, start, end):
		self._fid = fid
		self._chr = chr
		self._start = start
		self._end = end
		self._len = abs(start - end)

	@property
	def chr(self):
		return self._chr

	@property
	def start(self):
		return self._start

	@property
	def end(self):
		return self._end

	@property
	def len(self):
		return self._len


def compare_true_vs_reconstruct_single(ftrue,freconstruct):
	"""
	True vs reconstruct
	"""
	pd1 = pd.read_csv(ftrue, header=0, sep="\t")
	pd2 = pd.read_csv(freconstruct, header=0, sep="\t")

	pass


def similarity(f1, f2, props, match, penalty):
	# same fragment
	if f1 == f2:
		# length
		return match

	# check degree of overlap between fragments
	if props[f1].chr == props[f2].chr:
		if props[f2].start <= props[f1].start <= props[f2].end:
			# min (start_c1-end_c1,start_c1-end_c2)
			return max(abs(props[f1].start - props[f1].end), abs(props[f1].start - props[f2].end)) / props[f1].len
		if props[f1].start <= props[f2].start <= props[f1].end:
			# max (start_c2-end_c2,start_c2-end_c1)
			return max(abs(props[f2].start - props[f2].end), abs(props[f2].start - props[f1].end)) / props[f1].len

	return -2 * penalty


def score(m, i, j, i_props, j_props, sc):

	if i == 0 and j == 0:
		return sc

	if i == 0:
		return score(m, i, j - 1, i_props, j_props, sc + 0)

	if j == 0:
		return score(m, i - 1, j, i_props, j_props, sc + 0)

	# match
	if m[i - 1, j - 1] >= m[i - 1, j] and m[i - 1, j - 1] >= m[i, j - 1]:
		if m[i, j] >= 0:
			return score(m, i - 1, j - 1, i_props, j_props, sc + 1 * min(i_props[i - 1], j_props[j - 1]))
		else:
			return score(m, i - 1, j - 1, i_props, j_props, sc + 0)

	if m[i - 1, j] > m[i - 1, j - 1] and m[i - 1, j] > m[i, j - 1]:
		return score(m, i - 1, j, i_props, j_props, sc + 0)

	return score(m, i, j - 1, i_props, j_props, sc + 0)


def nw(x, y, props, match=1, penalty=1):
	nx = len(x)
	ny = len(y)

	if nx == 1 and ny == 1 and abs(x[0]) != abs(y[0]):
		return 0

	if nx == 1 and ny > 1:
		return 0

	if ny == 1 and nx > 1:
		return 0

	x_props = np.array([props[abs(f)].len for f in x])
	x_len = np.sum(x_props)
	x_props = x_props / x_len

	y_props = np.array([props[abs(f)].len for f in y])
	y_len = np.sum(y_props)
	y_props = y_props / y_len

	# Initialization
	m = np.zeros((nx + 1, ny + 1))
	m[0][0] = -1
	# rows
	for i in range(1, nx + 1):
		m[i, 0] = m[i - 1, 0] - 1
	# cols
	for i in range(1, ny + 1):
		m[0, i] = m[0, i - 1] - 1

	# Fill matrix
	for i in range(1, nx + 1):
		for j in range(1, ny + 1):
			m[i, j] = max(m[i - 1, j] - penalty,
			              m[i, j - 1] - penalty,
			              m[i - 1, j - 1] + similarity(abs(x[i - 1]), abs(y[j - 1]), props, match, penalty))

	# Compute score
	s = score(m, nx, ny, x_props, y_props, 0)
	return s

def compute_jaccard(c1, c2, graph):
	"""
	Compute Jaccard Index for 2 cycles in the graph and weight by the length
	Ignore orientation of the fragments

	Arguments:

	"""

	f1 = defaultdict(int)
	f2 = defaultdict(int)
	s1 = 0
	s2 = 0

	# count how many times a fragment is present in the path
	for f in c1:
		f1[abs(f)] += 1

	for f in c2:
		f2[abs(f)] += 1

	unique_fragments = list(set(list(f1.keys()) + list(f2.keys())))

	n_total = 0
	n_intersect = 0

	# check the overlap between circle c1 and c2
	for fid in unique_fragments:
		n_total += max(f1[fid], f2[fid]) * graph.fragments[fid].len
		n_intersect += min(f1[fid], f2[fid]) * graph.fragments[fid].len
		s1 += f1[fid] * graph.fragments[fid].len
		s2 += f2[fid] * graph.fragments[fid].len

	# jaccard_index, length_cycle1, length_cycle2
	return n_intersect / n_total, s1, s2

def compute_jaccard_unique_fragments(c1, c2, graph):
	"""
	Compute Jaccard Index for 2 cycles in the graph and weight by the length
	Ignore orientation of the fragments

	Arguments:

	"""

	f1 = defaultdict(int)
	f2 = defaultdict(int)
	s1 = 0
	s2 = 0

	# count how many times a fragment is present in the path
	for f in c1:
		f1[abs(f)] += 1

	for f in c2:
		f2[abs(f)] += 1

	unique_fragments = list(set(list(f1.keys()) + list(f2.keys())))

	n_total = 0
	n_intersect = 0

	# check the overlap between circle c1 and c2
	for fid in unique_fragments:
		n_total += graph.fragments[fid].len
		n_intersect += graph.fragments[fid].len if f1[fid] >= 1 and f2[fid] >= 1 else 0
		s1 += graph.fragments[fid].len
		s2 += graph.fragments[fid].len

	# jaccard_index, length_cycle1, length_cycle2
	return n_intersect / n_total, s1, s2

def filter_out_similar_cycles(cycles_list, graph, similarity_threshold = 0.9):
	"""
	# Do not keep very similar cycles for the inference
	# Remove high identity cycles (keep only one)
	"""
	n_cycles = len(cycles_list)
	to_remove = np.zeros(n_cycles)
	for i in range(0,n_cycles):
		if to_remove[i] == 0:
			for j in range(i,n_cycles):
				if i != j and to_remove[j] == 0:
					# s = nw(cycles_list[i], cycles_list[j], graph.fragments, match=1, penalty=1)
					# jscore, s1, s2 = compute_jaccard(cycles_list[i], cycles_list[j], graph)
					# print(s, jscore, cycles_list[i], cycles_list[j])
					# if s > similarity_threshold or (s > similarity_threshold / 2 and jscore > similarity_threshold):
					jscore, s1, s2 = compute_jaccard_unique_fragments(cycles_list[i], cycles_list[j], graph)
					if jscore > similarity_threshold:
						if s1 > s2:
							to_remove[j] = 1
						else:
							to_remove[i] = 1

	return [c for i,c in enumerate(cycles_list) if to_remove[i] == 0]


def bedtocycles_reconstructed(bedfile):
	"""
	Transform bed conformation to cycle list

	Arguments:
		bedfile (str): File

	Returns:
		cycles {c1: [(chr, start, stop, strand, frag, circid, cov),]},
		cycles_id [c1,c2,c3,...]
	"""
	cycles = {}
	df = pd.read_csv(bedfile, header=0, sep="\t")
	circ_ids = df[bp.CIRC_ID].drop_duplicates().tolist()

	for c in circ_ids:
		df_slice = df[df[bp.CIRC_ID] == c]
		cycles[c] = {}
		cycles[c]["cycles"] = []
		for i in range(0,df_slice.shape[0]):
			#chr	start	end	circ_id	fragment_id	strand	coverage	score
			row = df_slice.iloc[i]
			cycles[c]["cycles"].append((row[0],
			                  row[1],
			                  row[2],
			                  row[5],
			                  row[4],
			                  row[3],
			                  row[6]
			                  ))
	return cycles, circ_ids

def bedtocycles_true(bedfile):
	"""
	Transform bed conformation to cycle list

	Arguments:
		bedfile (str): File

	Returns:
		cycles {c1: [(chr, start, stop, strand, frag, circid, cov),]},
		cycles_id [c1,c2,c3,...]
	"""
	cycles = {}
	df = pd.read_csv(bedfile, header=0, sep="\t")
	circ_ids = df["structure"].drop_duplicates().tolist()

	for c in circ_ids:
		df_slice = df[df["structure"] == c]
		cycles[c] = {}
		cycles[c]["cycles"] = []
		for i in range(0,df_slice.shape[0]):
			#chr	start	stop	direction	target	coverage	structure	fragment
			row = df_slice.iloc[i]
			cycles[c]["cycles"].append((row[0],
			                  row[1],
			                  row[2],
			                  row[3],
			                  row[7],
			                  row[6],
			                  row[5]
			                  ))
	return cycles, circ_ids


def pident_regions(df, leng):
	"""
	Pidentity regions...

	Arguments:
		df (pd.DataFrame): Dataframe containing start and stop of matches
		leng (int): Length of the entire region
	"""
	if df.shape[0] == 0:
		return 0

	intervals_ = list(df.itertuples(index=False, name=None))
	tree_matches = IntervalTree.from_tuples(intervals_)

	# slice entire region based on the breakpoints
	tree_complete = IntervalTree.from_tuples([(1, leng)])
	breakpoints_subject = list(set(df.iloc[:, 0].tolist() + df.iloc[:, 1].tolist()))

	for b in breakpoints_subject:
		tree_complete.slice(b)
	for b in breakpoints_subject:
		tree_matches.slice(b)

	umatch = sorted(tree_complete - tree_matches)

	ulen = 0
	if len(umatch) > 0:
		for i in umatch:
			ulen += abs(i.begin - i.end)

	return 1 - ulen / leng


def compute_similarity(fasta1, fasta2, len1, len2, pident_threshold=90, bitscore_threshold=1000):

	# map against fasta1 (true)
	out_blast = subprocess.Popen(
		"""blastn -query {} -subject {} -outfmt 6 -ungapped -task megablast""".format(fasta2, fasta1), shell=True,
		stdout=subprocess.PIPE)
	# out_blast = subprocess.Popen("blastn -query s2.fasta -subject s1.fasta -outfmt 6 -gapopen 2 -gapextend 1",shell=True, stdout=subprocess.PIPE)
	df = pd.read_csv(StringIO(out_blast.communicate()[0].decode('utf-8')), header=None, sep="\t")
	df.columns = ["#qseqid", "sseqid", "pident", "alignment_length", "number_mismatches", "number_gap_openings",
	              "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]

	# keep only highly similar sequence identity
	df = df[(df["pident"] >= pident_threshold) & (df["bitscore"] >= bitscore_threshold)]
	# filter out inverted matches (we aim to search for correct matches to compute identity)
	df = df[df["subject_start"]<df["subject_end"]]
	pident1 = pident_regions(df[["subject_start", "subject_end"]], len1)

	# map against fasta2 (reconstruct)
	out_blast = subprocess.Popen(
		"""blastn -query {} -subject {} -outfmt 6 -ungapped -task megablast""".format(fasta1, fasta2), shell=True,
		stdout=subprocess.PIPE)
	# out_blast = subprocess.Popen("blastn -query s2.fasta -subject s1.fasta -outfmt 6 -gapopen 2 -gapextend 1",shell=True, stdout=subprocess.PIPE)
	df = pd.read_csv(StringIO(out_blast.communicate()[0].decode('utf-8')), header=None, sep="\t")
	df.columns = ["#qseqid", "sseqid", "pident", "alignment_length", "number_mismatches", "number_gap_openings",
	              "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]

	df = df[(df["pident"] >= pident_threshold) & (df["bitscore"] >= bitscore_threshold)]
	df = df[df["subject_start"] < df["subject_end"]]
	pident2 = pident_regions(df[["subject_start", "subject_end"]], len2)

	return pident1, pident2


def get_fasta_perinvertal(_chr, _start, _stop, _strand, circ_id, ref):
	a = pybedtools.BedTool("""{} {} {} {} . {}""".format(str(_chr),
	                                                     str(_start),
	                                                     str(_stop),
	                                                     str(circ_id),
	                                                     _strand), from_string=True)
	a = a.sequence(fi=ref, s=True)
	return open(a.seqfn).read().split('\n')[1]


def get_fasta(intervals, circ_id, ref):
	"""
	Get fasta based on the bed regions
	"""
	s = ""
	for i in intervals:
		_chr = i[0]
		_start = i[1]
		_stop = i[2]
		_strand = i[3]
		s += get_fasta_perinvertal(_chr, _start, _stop, _strand, circ_id, ref)

	return s

def compare_true_reconstruct(c_true, c_reconstruct, ref, out, temp="/temp"):
	"""
	Compare true circular structure against reconstructed

	Arguments:
		c_true (bed): Configuration file
		c_reconstruct (bed): Configuration file reconstruction
		temp (str): Temp root
		outfile (txt): Outputfile
		           true_id true_len true_topology reconstruct_id reconstruct_len reconstruct_topology pident_true pident_reconstruct
	"""
	true_cycles, true_circ_id_list = bedtocycles_true(c_true)
	recons_cycles, recons_circ_id_list = bedtocycles_reconstructed(c_reconstruct)

	for c in true_circ_id_list:
		true_cycles[c]["topology"], true_cycles[c]["topology_details"] = metrics.annotate_topology(true_cycles[c]["cycles"])

	for c in recons_circ_id_list:
		recons_cycles[c]["topology"], recons_cycles[c]["topology_details"] = metrics.annotate_topology(recons_cycles[c]["cycles"])

	pprint.pprint(true_cycles)
	pprint.pprint(recons_cycles)

	with open(out, "w") as f:

		f.write("true_id\ttrue_len\ttrue_topology\treconstruct_id\treconstruct_len\treconstruct_topology\tpident_true\tpident_reconstruct\n")
		for c1 in true_circ_id_list:
			for c2 in recons_circ_id_list:

				s1 = get_fasta(true_cycles[c1]["cycles"], c1, ref)
				s2 = get_fasta(recons_cycles[c2]["cycles"], c2, ref)

				rec1 = SeqRecord(Seq(s1), id=str(c1), name=str(c1))
				rec2 = SeqRecord(Seq(s2), id=str(c2), name=str(c2))

				os.makedirs(temp,exist_ok=True)
				p1 = os.path.join(temp,str(uuid.uuid4()) + "_s1.fasta")
				p2 = os.path.join(temp,str(uuid.uuid4()) + "_s2.fasta")
				SeqIO.write([rec1], p1, "fasta")
				SeqIO.write([rec2], p2, "fasta")

				(pident1, pident2) = compute_similarity(p1, p2, len(s1), len(s2))
				f.write("""{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n""".format(str(c1),
				                                      str(len(s1)),
				                                      str(true_cycles[c1]["topology"]),
				                                      str(c2),
				                                      str(len(s2)),
				                                      str(recons_cycles[c2]["topology"]),
				                                      str(pident1),
				                                      str(pident2)))

def transform_rawsim_to_topology_details(jdata):
	"""
	Transform the input array for a simulated ecDNA to topology_details format
	"""
	jnew = {}
	header = ["#frag", "len", "small_del", "dup", "inv", "interch", "multi+region", "foldback", "return"]
	for i in range(0,len(jdata["accept_criteria"])):
		jnew[header[i]] = jdata["accept_criteria"][i]
	
	return jnew





