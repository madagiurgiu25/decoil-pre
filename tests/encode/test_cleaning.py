import unittest
import os
import re
import pysam
import logging
import numpy as np

import decoil.encode.operations as operations
import decoil.encode.encode as encode
import decoil.search.cycles as cycles
from decoil.encode.graph import MultiGraph
from decoil.utils import VCF_PROP as vp

logger = logging.getLogger("decoil.encode")


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
  

def test_clean_multi_vcf(self):
	""" Test multi-VCF file
	"""

def test_read_and_clean_vcf(self):
	os.makedirs("tests/examples/vcfs/output_sim1010_sniffles1_withartifacts", exist_ok=True)
	encode.readvcf("tests/examples/vcfs/sim1010_sniffles1_withartifacts.vcf", 
					"tests/examples/vcfs/output_sim1010_sniffles1_withartifacts", multi=False, svcaller=vp.SNIFFLES1)
	# test if clean.vcf and clean_filtered.vcf exists
	
	# test if the clean.vcf the same with input
	
	# test if the clean_filtered.vcf removed specific columns

# pytest injects it
def test_clean_corrupted_vcf_warning(caplog):
	with caplog.at_level(logging.DEBUG, logger="decoil.encode"):
		os.makedirs("tests/examples/vcfs/output_sim1010_sniffles1_corrupted", exist_ok=True)
		encode.cleanvcf("tests/examples/vcfs/sim1010_sniffles1_corrupted.vcf", 
						os.path.join("tests/examples/vcfs/output_sim1010_sniffles1_corrupted","clean.vcf"),
						multi=False)
	warning_occurrences = re.findall(r"WARNING", str(caplog.text), flags=re.MULTILINE)
	count_warnings = len(warning_occurrences)
	assert 2 == count_warnings

if __name__ == '__main__':
	unittest.main()
