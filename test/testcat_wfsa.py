import cat_wfsa
from math import log
from wfsa import integer_code_len
import unittest

class TestCatWFSA(unittest.TestCase):
	def setUp(self):
		self.alphabet = set('aeiouptkbdgszhlmn')
		self.categories = \
		{
			'V':set('aeiou'),
			'C':set('ptkbdgszhlmn'),
			'Son':set('lmn'),
			'Nas':set('mn'),
			'Obs':set('ptkbdgszh'),
			'Stop':set('ptkbdg'),
			'Fric':set('szh'),
			'Cont':set('szhlmn'),
			'Voice':set('bdgzlmn'),
			'Vless':set('ptksh')
		}
		self.word = 'abstlakt'
		self.voice_voiceless = cat_wfsa.CatWFSA(
			'aeiouptkbdgszhlmn', '#', {'#':-1,'VStop':-1}, 
			{
				'#':{frozenset(['Voice','Stop']):('VStop', -1),'_other':('#', -1)},
				'VStop':{frozenset(['Vless','Fric']):('#', 0),'_other':('#',-1)}
			},
			self.categories, [1.5], False )
		
		self.ccvcc = cat_wfsa.CatWFSA(
			'aeiouptkbdgszhlmn', '#', dict([(x,-1) for x in '#1234']),
			{
				'#':{frozenset(['C']):('1', -1),'_other':('#', -1)},
				'1':{frozenset(['C']):('2', -1),'_other':('#', -1)},
				'2':{frozenset(['V']):('3', -1),frozenset(['C']):('2', -1)},
				'3':{frozenset(['C']):('4', -1),'_other':('#', -1)},
				'4':{frozenset(['C']):('2', 0),'_other':('#', 0)}
			},
			self.categories, [3], False )
	
	def test_weight(self):
		self.assertAlmostEqual(self.voice_voiceless.weight(self.word), 1.5)
		self.assertAlmostEqual(self.voice_voiceless.weight('slap'), 0.0)
		self.assertAlmostEqual(self.ccvcc.weight(self.word), 3)
		self.assertAlmostEqual(self.ccvcc.weight('ppeskutt'), 6)
	
	def test_intersect(self):
		intersect = self.voice_voiceless.intersect(self.ccvcc)
		self.assertAlmostEqual(intersect.weight(self.word), 4.5)
	
	def test_complexity(self):
		#2 states, 18 letters, 10 categories, 4 arcs, 2 with 2 categorios,
		#2 with a letter. 2 stop states, 1 parameter
		two = integer_code_len(2)
		arcs = integer_code_len(4)
		arc_labels = 2*(two+2*log(10,2)) + 2*log(18, 2)
		arc_encodings = arcs + 12 + arc_labels
		complexity = 2*two + arc_encodings + 4 + integer_code_len(1) + 64
		self.assertAlmostEqual(complexity, self.voice_voiceless.complexity)
	
if __name__ == "__main__":
	unittest.main()