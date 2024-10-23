"""
Created using PyCharm
@author: Madalina Giurgiu
@date: 12:11 PM 9/20/22

Prepare reconstruction output.
"""
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

from decoil.utils import BED_PROP as bp
from decoil.utils import LINKS_PROP as lp
from decoil.validate import compare as compare
from decoil.utils import TOPOLOGY as tp

def to_fragments_list(df, cid, counter):
	"""
	Input header: #chr	start	end	circ_id	fragment_id	strand	coverage	score

	Returns fragment list
	"""
	df_slice = df[df.circ_id == cid]
	frags = {}
	frags_list = []
	for i in range(df_slice.shape[0]):
		counter += 1
		id = counter
		_chr = df_slice.iloc[i,0]
		_start = df_slice.iloc[i, 1]
		_stop = df_slice.iloc[i, 2]
		frags[id] = compare.Prop(id,_chr,_start,_stop)
		frags_list.append(id)

	return frags_list, frags, counter


def convert_bed2links(bedfile, linksfile):
	"""
	For Ggnome visualization we need a specific format for the links between connected fragments
	"""
	dict_links = []  # array of dicts
	
	df = pd.read_csv(bedfile, header=0, sep="\t")
	circles = df[bp.CIRC_ID].drop_duplicates().tolist()
	
	for circ in circles:
		df_temp = df.loc[df[bp.CIRC_ID] == circ]
		lencirc = df_temp.shape[0]
		
		# print(df_temp)
		# print(lencirc)
		
		for i in range(0, lencirc):
			# print(i)
			pos1 = i % lencirc
			pos2 = (i + 1) % lencirc
			
			n1 = df_temp.index[pos1]
			n2 = df_temp.index[pos2]
			
			strand1 = df_temp.iloc[pos1, bp.STRAND_IDX]
			strand2 = df_temp.iloc[pos2, bp.STRAND_IDX]
			
			# tail/head equivalent with left/right in this annotation
			if strand1 == bp.POSITIVE and strand2 == bp.POSITIVE:
				# connect head/right to tail/left
				n1side = bp.RIGHT
				n2side = bp.LEFT
			elif strand1 == bp.NEGATIVE and strand2 == bp.POSITIVE:
				# connect tail/left to tail/left
				n1side = bp.LEFT
				n2side = bp.LEFT
			elif strand1 == bp.POSITIVE and strand2 == bp.NEGATIVE:
				# connect head/right to head/right
				n1side = bp.RIGHT
				n2side = bp.RIGHT
			else:
				# connect tail/left to head/right
				n1side = bp.LEFT
				n2side = bp.RIGHT
			
			elem = {lp.N1: n1 + 1,  # 1-based instead of 0-based
			        lp.N2: n2 + 1,
			        lp.N1SIDE: n1side,
			        lp.N2SIDE: n2side,
			        lp.CIRC_ID: circ}
			dict_links.append(elem)
	
	# save
	dfnew = pd.DataFrame(dict_links)
	dfnew.to_csv(linksfile, sep="\t", header=True, index=False)


def convert_cycles2bed(paths, bedfile, graph):
	"""
	Convert cycles to bed
	"""
	fragments = graph.get_fragments()

	with open(bedfile, "w") as f:
		f.write("#chr\tstart\tend\tcirc_id\tfragment_id\tstrand\tcoverage\testimated_proportions\n")
		for i, path in enumerate(paths):
			for fid in path:
				# fid can be negative denoting reverse strand
				absfid = abs(fid)
				chr, start, stop = fragments[absfid].chr, fragments[absfid].start, fragments[absfid].end
				circid = "circ_{}".format(str(i))
				frag = absfid
				strand = "+" if fid > 0 else "-"
				cov = fragments[absfid].coverage
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(chr),
				                                              str(start),
				                                              str(stop),
				                                              str(circid),
				                                              str(frag),
				                                              strand,
				                                              str(int(cov)),
				                                              "0"))



def convert_path2bed(paths, bedfile):
	"""
	Convert graph paths to bed format
	"""
	
	with open(bedfile, "w") as f:
		f.write("#chr\tstart\tend\tcirc_id\tfragment_id\tstrand\tcoverage\testimated_proportions\n")
		for path in paths:
			for (chr, start, stop, strand, frag, circid, cov) in path:
				f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(chr),
				                                              str(start),
				                                              str(stop),
				                                              str(circid),
				                                              str(frag),
				                                              strand,
				                                              str(int(cov)),
				                                              "0"))


def convert_path2bed_v2(path, bedfile, keep=None, score_threshold=0):
	"""
	Convert graph path to bed format.

	Arguments:
		path (dict): Dict storing the circular paths candidates.
					 Example:
						{7: {'conf': [('1', '+', 'a:h:a:t:0', 50),
									 ('2', '+', 'b:h:b:t:0', 20)],
						  'score': 35.0},
						 4: {'conf': [('3', '+', 'c:h:c:t:0', 110),
									   ('1', '+', 'a:h:a:t:0', 50),
									   ('2', '+', 'b:h:b:t:0', 20),
						   'score': 28.78723617397549}
						}
		bedfile (str):
		keep (str): Keep ecDNA annotated circles or all
		score_threshold (int): Keep only circles with a minimal score (estimated copy-number)
	"""
	df = pd.DataFrame(columns=["#chr","start","end","circ_id","fragment_id","strand","coverage","estimated_proportions"])

	to_keep_circles = {}
	for id_circ in path:
		# skip non ecDNA entries
		if keep == "ecDNA" and path[id_circ]["label"] != "ecDNA":
			continue
		score = path[id_circ]["score"]
		if int(score) >= score_threshold:
			to_keep_circles[id_circ] = int(score)

	# sort circles based on the estimated proportions (decreasing
	sorted_dict = dict(sorted(to_keep_circles.items(), key=lambda item: -item[1]))

	for id_circ in sorted_dict:
		for (chr, start, stop, strand, frag, circid, cov) in path[id_circ]["conf"]:
			new_row = {'#chr': str(chr),
					   'start': str(start),
					   'end': str(stop),
					   'circ_id': str(circid),
					   'fragment_id': str(frag),
					   'strand': strand,
					   'coverage': str(int(cov)),
					   'estimated_proportions': str(int(path[id_circ]["score"]))
					   }
			df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

	df.to_csv(bedfile,sep="\t",header=True,index=False,quoting=None)

def convert_bed2fasta(bedfile, fastafile, ref_genome, version=1):
	"""
	Convert a graph path bed to fasta.
	"""
	
	# parse faste file and turn into dictionary
	records = SeqIO.to_dict(SeqIO.parse(open(ref_genome), 'fasta'))
	
	circleseq_records = []
	sequence_circle = ""
	lastid = -1
	
	with open(bedfile) as f:
		for line in f:
			if line.startswith("#"):
				continue
			
			# different computation versions
			if version == 1:
				chr, start, stop, strand, fragnr, circleid, cov = line.split()
			else:
				chr, start, stop, circleid, fragnr, strand, cov, score = line.split()
			
			start = int(start)
			stop = int(stop)
			
			if lastid == -1:
				lastid = circleid
			
			# if a complete circle was constructed
			if lastid != circleid:
				circle_record = SeqRecord(Seq(sequence_circle), id=str(lastid), description=str(lastid))
				circleseq_records.append(circle_record)

				sequence_circle = ""
				lastid = circleid
			
			# get sequence
			long_seq_record = records[chr]
			long_seq = long_seq_record.seq
			
			short_seq = str(long_seq)[start:stop]
			if strand == "-":
				short_seq = str(Seq(short_seq).reverse_complement())
			sequence_circle += short_seq
	
	# last element
	circle_record = SeqRecord(Seq(sequence_circle), id=str(lastid), description=str(lastid))
	circleseq_records.append(circle_record)
	
	# write to file
	with open(fastafile, 'w') as f:
		SeqIO.write(circleseq_records, f, 'fasta')


def create_summary(candidates, summaryfile):
	df = pd.DataFrame(
		columns=["circ_id", "chr_origin", "size(MB)", "label", "topology_idx", "topology_name", "estimated_proportions"])

	for c in candidates:
		top_idx = candidates[c]["topology"]
		top_name = tp.DICT[top_idx]

		new_row = { "circ_id": str(c),
					"chr_origin": str(candidates[c]["chrs"]),
					"size(MB)": str(candidates[c]["size"]),
					"label": str(candidates[c]["label"]),
					"topology_idx": str(top_idx),
					"topology_name": top_name,
					"estimated_proportions": str(int(candidates[c]["score"]))
				   }
		# df = df.append(new_row, ignore_index=True)
		df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

	df.to_csv(summaryfile, sep="\t", header=True, index=False, quoting=None)

	#
	# with open(summaryfile, "w") as f:
	# 	f.write("circ_id\tchr_origin\tsize(MB)\tlabel\ttopology_idx\ttopology_name\testimated_proportions\n")
	# 	for c in candidates:
	# 		top_idx = candidates[c]["topology"]
	# 		top_name = tp.DICT[top_idx]
	# 		f.write("""{}\t{}\t{}\t{}\t{}\t{}\t{}\n""".format(str(c),
	# 		                                str(candidates[c]["chrs"]),
	# 		                                str(candidates[c]["size"]),
	# 		                                str(candidates[c]["label"]),
	# 		                                str(top_idx),
	# 		                                top_name,
	# 										str(int(candidates[c]["score"]))
	# 		                                ))



