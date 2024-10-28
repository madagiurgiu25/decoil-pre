import unittest

from decoil.encode.graph import MultiGraph
from decoil.utils import GRAPH_PROP as gp
import decoil.search.search as search

def cf(paths):
	d = {}
	d[1] = {2: (0, "A")}
	d[2] = {1: (0, "-A")}

	d[3] = {4: (2, "B")}
	d[4] = {3: (2, "-B")}

	d[5] = {6: (4, "C")}
	d[6] = {5: (4, "-C")}

	d[7] = {8: (6, "F")}
	d[8] = {7: (6, "-F")}

	for p in paths:
		s = ""
		for i in range(0, len(p) - 1):
			if p[i][0] in d and p[i + 1][0] in d[p[i][0]]:
				s += d[p[i][0]][p[i + 1][0]][1]
		print(s)

class DFSTest(unittest.TestCase):

	def start_search(self, graph, find_function):

		graph.reset_time()

		for comp in graph.connected_components():
			gcomp = graph.subgraph(comp)
			findings = []

			# 1. find circular path starting from each node of the component
			nodes = gcomp.get_nodes()
			for v in nodes:
				if nodes[v].part == gp.TAIL:
					print("From", nodes[v].parent_fragment)
					found_paths = find_function(v)
					cf(found_paths)
					print()
					findings = search.deduplicate(findings, found_paths)


	# def test_ABF_BCF(self):

	# 	# load the ABF, BCF graph
	# 	graph = MultiGraph.load_graph("tests/data/test2.txt")
	# 	self.assertEqual(len(graph.edges), 10)
	# 	self.assertEqual(graph.get_edges()[7].u, 3)
	# 	self.assertEqual(graph.get_edges()[7].v, 8)
	# 	self.assertEqual(graph.get_edges()[7].svtype, 2)

	# 	print("test_ABF_BCF_v1")
	# 	self.start_search(graph, graph.find_simple_circles)
	# 	print("test_ABF_BCF_v2")
	# 	self.start_search(graph, graph.find_simple_circles2)

if __name__ == '__main__':
	unittest.main()
