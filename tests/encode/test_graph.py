import unittest

import pysam
import numpy as np

import decoil.encode.operations as operations
import decoil.search.cycles as cycles
from decoil.encode.graph import MultiGraph


class TestMultiGraph(unittest.TestCase):

	# def test_adjancency_matrix(self):
	# 	g = MultiGraph.load_graph("tests/data/test1.txt")
	# 	m = g.adjacency_matrix()
	# 	m_compare = [[0, 1, 0, 1, 0, 0],
	# 	             [1, 0, 1, 0, 0, 0],
	# 	             [0, 0, 0, 1, 0, 1],
	# 	             [1, 0, 1, 0, 1, 0],
	# 	             [0, 0, 0, 0, 0, 1],
	# 	             [0, 0, 1, 0, 1, 0]
	# 	             ]
	# 	self.assertTrue(np.array_equal(m, m_compare, equal_nan=True))

	def test_insert_beginning_list(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 150,
		        "rank2": 50,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 7)
	
	def test_insert_middle_list(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 90,
		        "rank2": 50,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[1][0], 7)
	
	def test_insert_end_list(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 70,
		        "rank2": 50,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[2][0], 7)
	
	def test_insert_beginning_rank2_list(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 100,
		        "rank2": 60,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 7)
		self.assertEqual(newlist[1][0], 1)
		self.assertEqual(newlist[2][0], 2)
	
	def test_insert_beginning_rank2_list2(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 100,
		        "rank2": 40,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 1)
		self.assertEqual(newlist[1][0], 2)
		self.assertEqual(newlist[2][0], 7)
	
	def test_insert_beginning_rank2_list3(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 80,
		        "rank2": 60,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 1)
		self.assertEqual(newlist[1][0], 7)
		self.assertEqual(newlist[2][0], 2)
	
	def test_insert_beginning_rank2_list4(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 80,
		        "rank2": 40,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 1)
		self.assertEqual(newlist[1][0], 2)
		self.assertEqual(newlist[2][0], 7)
	
	def test_insert_beginning_rank3_list(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 100,
		        "rank2": 50,
		        "rank3": 2}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 7)
		self.assertEqual(newlist[1][0], 1)
		self.assertEqual(newlist[2][0], 2)
	
	def test_insert_beginning_rank3_list2(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 100,
		        "rank2": 50,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 1)
		self.assertEqual(newlist[1][0], 7)
		self.assertEqual(newlist[2][0], 2)
	
	def test_insert_beginning_rank3_list3(self):
		
		# test insert beginning list
		list_edges = [(1, {"rank1": 100,
		                   "rank2": 50,
		                   "rank3": 3}),
		              (2, {"rank1": 80,
		                   "rank2": 50,
		                   "rank3": 3})]
		e = 7
		data = {"rank1": 80,
		        "rank2": 50,
		        "rank3": 3}
		
		newlist = cycles.insert_edges_sorted(list_edges, e, data)
		self.assertEqual(newlist[0][0], 1)
		self.assertEqual(newlist[1][0], 2)
		self.assertEqual(newlist[2][0], 7)
	
	# def test_load_graph(self):
	# 	g = MultiGraph.load_graph("tests/data/test1.txt")
	# 	self.assertEqual(len(g.edges), 9)
	# 	self.assertEqual(g.get_edges()[8].u, 3)
	# 	self.assertEqual(g.get_edges()[8].v, 6)
	# 	self.assertEqual(g.get_edges()[8].svtype, 1)
	
	def test_merge_breakpoints(self):
		
		breakpoints = {"chrM": [0, 16550, 16570]}
		svinfo = {"chrM_0": [('60', 'chrM', '16550', 'DUP', 9, 0, '1/1', '-+'),
		                     ('50', 'chrM', '16570', 'DUP', 9, 0, '1/1', '-+')]}
		collection_breakpoints, svinfo = operations.merge_near_breakpoints(breakpoints, svinfo)
		print(collection_breakpoints, svinfo)
		self.assertEqual(len(breakpoints["chrM"]), 2)
	
	def test_merge_breakpoints2(self):
		
		breakpoints = {"chrM": [0, 16550, 16570, 16680]}
		svinfo = {"chrM_0": [('60', 'chrM', '16550', 'DUP', 9, 0, '1/1', '-+'),
		                     ('50', 'chrM', '16570', 'DUP', 9, 0, '1/1', '-+')],
		          "chrM_16570": [('50', 'chrM', '16680', 'DEL', 9, 0, '1/1', '-+')]}
		
		collection_breakpoints, svinfo = operations.merge_near_breakpoints(breakpoints, svinfo)
		print(collection_breakpoints, svinfo)
	
	# def test_parse_read(self):
	# 	samfile = pysam.AlignmentFile("tests/data/ngmlr.bam", "rb")
	# 	breakpoint = 15102001
	# 	window = 500
	# 	padd = 100
	# 	span = 500000
	# 	count_coverage_through = 0
		
	# 	read_aligned = ["S1_1743", "S1_886"]
	# 	found = 0
	# 	read_split = ["S1_1389"]
		
	# 	count_total = samfile.count('chr2', breakpoint - window, breakpoint + window)
		
	# 	for read in samfile.fetch(contig='chr2',
	# 	                          start=breakpoint - window,
	# 	                          stop=breakpoint + window):
			
	# 		# filter out secondary alignment and mapped
	# 		if not read.is_mapped or read.is_secondary:
	# 			continue
			
	# 		# filter out short reads (< 2*window)
	# 		if read.query_length < 2 * window:
	# 			continue
			
	# 		# keep only reads which go through the breakpoint
	# 		if read.reference_start < (breakpoint - padd) and read.reference_end > (breakpoint + padd):
	# 			# print("through read",read.query_name,read.query_alignment_start, read.reference_start)
	# 			if read.query_name in read_aligned:
	# 				found += 1
	# 			# print("Found in read_aligned")
	# 			count_coverage_through += 1
		
	# 	self.assertEqual(found, 2)
	# 	self.assertEqual(count_coverage_through, 40)
	# 	self.assertEqual(count_total, 98)

# 	count_coverage_through += 1
# 	# check if the breakpoint split the read
# 	total = read.query_length
# 	start = read.reference_start
# 	stop = read.reference_end
# 	reflen = read.query_alignment_length
# 	alstart = read.query_alignment_start
# 	alstop = read.query_alignment_end
# 	split = False
# 	# if (breakpoint - window) > start or (breakpoint + window) > start
# 	print(read.query_name, alstart, alstop, total, (alstart+alstop)/total, abs(alstart-alstop), reflen, read.is_forward)
#

# print("Coverage", count_coverage_through)


if __name__ == '__main__':
	unittest.main()
