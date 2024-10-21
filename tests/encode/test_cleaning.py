import unittest

import pysam
import numpy as np

import decoil.encode.operations as operations
import decoil.search.cycles as cycles
from decoil.encode.graph import MultiGraph


class TestMultiGraph(unittest.TestCase):

	def test_remove_lowcoverage(self):
		g = MultiGraph.load_graph("tests/data/test3.txt")
		m = g.adjacency_matrix()
		MultiGraph.id_generator = 100
		m_initial = np.array([[0, 1, 0, 0, 0, 0, 0, 1],
		             		  [1, 0, 1, 0, 1, 0, 0, 0],
		                      [0, 1, 0, 1, 0, 0, 0, 1],
		                      [0, 0, 1, 0, 1, 0, 0, 0],
		                      [0, 1, 0, 1, 0, 1, 0, 0],
		                      [0, 0, 0, 0, 1, 0, 1, 0],
		                      [0, 0, 0, 0, 0, 1, 0, 1],
		                      [1, 0, 1, 0, 0, 0, 1, 0]
		             ])
		self.assertTrue(m.dtype, m_initial.dtype)
		self.assertTrue(np.array_equal(m, m_initial, equal_nan=True))

		# remove low coverage fragments
		g = operations.remove_lowcoverage_fragments(g, threshold=6)
		m = g.adjacency_matrix()
		fragments = g.get_fragments().keys()
		m_aftercleaning = np.array([[0, 1, 0, 0, 0, 1],
		             	   			[1, 0, 1, 0, 0, 0],
		                   			[0, 1, 0, 1, 0, 0],
		                   			[0, 0, 1, 0, 1, 0],
		                   			[0, 0, 0, 1, 0, 1],
		                   			[1, 0, 0, 0, 1, 0]])

		self.assertTrue(fragments, ['A','C','D'])
		self.assertTrue(np.array_equal(m, m_aftercleaning, equal_nan=True))

	def test_remove_shortfragments(self):
		g = MultiGraph.load_graph("tests/data/test3.txt")
		g.fragments['A'].len = 1000
		g.fragments['B'].len = 100
		g.fragments['C'].len = 20000
		g.fragments['D'].len = 1500

		m = g.adjacency_matrix()
		MultiGraph.id_generator = 100
		m_initial = [[0, 1, 0, 0, 0, 0, 0, 1],
		             [1, 0, 1, 0, 1, 0, 0, 0],
		             [0, 1, 0, 1, 0, 0, 0, 1],
		             [0, 0, 1, 0, 1, 0, 0, 0],
		             [0, 1, 0, 1, 0, 1, 0, 0],
		             [0, 0, 0, 0, 1, 0, 1, 0],
		             [0, 0, 0, 0, 0, 1, 0, 1],
		             [1, 0, 1, 0, 0, 0, 1, 0]
		             ]
		self.assertTrue(np.array_equal(m, m_initial, equal_nan=True))

		# remove low coverage fragments
		g = operations.remove_short_fragments(g, threshold=500)
		m = g.adjacency_matrix()
		fragments = g.get_fragments().keys()
		m_aftercleaning = [[0, 1, 0, 0, 0, 1],
		             [1, 0, 1, 0, 0, 0],
		             [0, 1, 0, 1, 0, 1],
		             [0, 0, 1, 0, 1, 0],
		             [0, 0, 0, 1, 0, 1],
		             [1, 0, 1, 0, 1, 0]
		             ]
		self.assertTrue(fragments, ['A','C','D'])
		self.assertTrue(np.array_equal(m, m_aftercleaning, equal_nan=True))

if __name__ == '__main__':
	unittest.main()
