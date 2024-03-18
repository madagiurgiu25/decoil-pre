import unittest

from decoil.search import cycles

class TestCyclesModule(unittest.TestCase):

	def test_rotate_to_minimal(self):
		"""
		Test if rotation works
		"""
		rotation = cycles.rotate_tominimal([3,2,4,5])
		self.assertEqual(rotation, [2,4,5,3])

	def test_rotate_to_minimal2(self):
		"""
		Test if rotation works
		"""
		rotation = cycles.rotate_tominimal([3,3,2,4,5])
		self.assertEqual(rotation, [2,4,5,3,3])

	def test_rotate_to_minimal3(self):
		"""
		Test if rotation works
		"""
		rotation = cycles.rotate_tominimal([3,3,2,2,4,5])
		self.assertEqual(rotation, [2,2,4,5,3,3])

	def test_rotate_to_minimal_negative(self):
		"""
		Test if rotation works
		"""
		rotation = cycles.rotate_tominimal([3,2,4,-5])
		self.assertEqual(rotation, [2,4,-5,3])

	def test_rotate_to_minimal_negative2(self):
		"""
		Test if rotation works
		"""
		rotation = cycles.rotate_tominimal([3,-2,2,4,-5])
		self.assertEqual(rotation, [-2,2,4,-5,3])

	def test_canonical_path1(self):
		"""
		Test canonical path conformation
		"""
		canonical = cycles.convert2canonical_path([4,3,2,1])
		self.assertEqual(canonical, [1,4,3,2])

	def test_canonical_path2(self):
		"""
		Test canonical path conformation
		"""
		canonical = cycles.convert2canonical_path([4,-3,-2,1])
		self.assertEqual(canonical, [1,4,-3,-2])

	def test_canonical_path3(self):
		"""
		Test canonical path conformation
		"""
		canonical = cycles.convert2canonical_path([4,-3,-2,-1,5])
		self.assertEqual(canonical, [1,2,3,-4,-5])

if __name__ == '__main__':
	unittest.main()
