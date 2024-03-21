import unittest

from decoil.encode.graph import MultiGraph

class TestGraphCycles(unittest.TestCase):

	def test_hash_match(self):
		m = MultiGraph()
		metahash1 = m.__metahash_circle__([1,2,3])
		metahash2 = m.__metahash_circle__([2,3,1])
		self.assertEqual(metahash1, metahash2)

	def test_hash_unmatch(self):
		m = MultiGraph()
		metahash3 = m.__metahash_circle__([2,3,1,4])
		metahash4 = m.__metahash_circle__([2,3,4,1])
		self.assertNotEqual(metahash3, metahash4)

	def test_hash_reverse_match(self):
		m = MultiGraph()
		metahash1 = m.__metahash_circle__([1,2,3])
		metahash2 = m.__metahash_circle__([-2,-1,-3])
		self.assertEqual(metahash1, metahash2)

	def test_hash_reverse_unmatch(self):
		m = MultiGraph()
		metahash1 = m.__metahash_circle__([1,2,3])
		metahash2 = m.__metahash_circle__([-2,-3,-1]) # this is not correct
		self.assertNotEqual(metahash1, metahash2)

if __name__ == '__main__':
	unittest.main()
