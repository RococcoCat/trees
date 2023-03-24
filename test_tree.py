import unittest
import trees
import pandas as pd
import numpy as np
import math
import sys

class TestSum(unittest.TestCase):

	def test_jukes_cantor(self):
		seq1 = "AAAAA"
		seq2 = "ACGTA"
		x = round(trees.jc_dist(seq1, seq2), 4)
		self.assertEqual(x, 1.2071)

	def test_check_ultrametric(self):
		df1 = pd.DataFrame(data = [[0, 1, 2], [1, 0, 3], [2, 3, 0]])
		df2 = pd.DataFrame(data = [[0, 3, 8, 8], [3, 0, 8, 8], [8, 8, 0, 3], [8, 8, 3, 0]])
		self.assertEqual(trees.is_ultrametric(df1), False)
		self.assertEqual(trees.is_ultrametric(df2), True)

	def test_check_upgma(self):
		cols = ['dog', 'bear', 'raccoon', 'weasel', 'seal', 'sea_lion', 'cat', 'monkey']
		dog = [0, 32, 48, 51, 50, 48, 98, 148]
		bear = [32, 0, 26, 34, 29, 33, 84, 136]
		raccoon = [48, 26, 0, 42, 44, 44, 92, 152]
		weasel = [51, 34, 42, 0, 44, 38, 86, 142]
		seal = [50, 29, 44, 44, 0, 24, 89, 142]
		sea_lion = [48, 33, 44, 38, 24, 0, 90, 142]
		cat = [98, 84, 92, 86, 89, 90, 0, 148]
		monkey = [148, 136, 152, 142, 142, 142, 148, 0]
		x = [dog, bear, raccoon, weasel, seal, sea_lion, cat, monkey]
		df = pd.DataFrame(data = x, columns = cols, index = cols)
		ans = '((dog:23.875, (((bear:13.0, raccoon:13.0):18.75, (seal:12.0, sea_lion:12.0):18.75):19.75, weasel:19.75):23.875):46.34375, cat:46.34375), monkey)'
		self.assertEqual(trees.upgma(df, "test.txt", True),ans)

	def test_check_neighbor_joining(self):
		pass

if __name__ == '__main__':
    unittest.main()