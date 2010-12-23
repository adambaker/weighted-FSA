from nat_class_wfsa import *
from nat_class_set import *
from math import log
from wfsa import integer_code_len
import unittest

class TestNatClassWFSA(unittest.TestCase):
	def setUp(self):
		self.alphabet = set('aeuptkbdszh')
		self.classes = set(
		[
			frozenset('aeu'),
			frozenset('ptkbdzsh'),
			frozenset('ptkbd'),
			frozenset('szh'),
			frozenset('bdz'),
			frozenset('ptksh'),
			frozenset('ptbdsz')
		])
		
		self.nat_class_set = NaturalClassSet( self.alphabet, self.classes )
		
		self.v_vless_states = [
			
		]
		
		self.voice_voiceless = NaturalClassWFSA(
			self.nat_class_set, '#', {'#':-1,'VStop':-1}, 
			{
				'#':{frozenset('bd'):('VStop', -1),'_other':('#', -1)},
				'VStop':{frozenset('sh'):('#', 0),'_other':('#',-1)}
			},
			[1.5], precision=False )
		
		self.ccvcc = NaturalClassWFSA(
			self.nat_class_set, '#', dict([(x,-1) for x in '#1234']),
			{
				'#':{frozenset('ptkbdszh'):('1', -1),'_other':('#', -1)},
				'1':{frozenset('ptkbdszh'):('2', -1),'_other':('#', -1)},
				'2':{frozenset('aeu'):('3', -1),frozenset('ptkbdszh'):('2', -1)},
				'3':{frozenset('ptkbdszh'):('4', -1),'_other':('#', -1)},
				'4':{frozenset('ptkbdszh'):('2', 0),'_other':('#', 0)}
			},
			[3], precision=False )

	def test_classes(self):
		classes = set(
		[
			frozenset('a'),frozenset('e'),frozenset('u'),frozenset('p'),
			frozenset('t'),frozenset('k'),frozenset('b'),frozenset('d'),
			frozenset('s'),frozenset('z'),frozenset('h'),
			frozenset('aeu'),frozenset('ptkbdszh'),frozenset('ptkbd'),
			frozenset('szh'),frozenset('ptksh'),frozenset('sh'),
			frozenset('bdz'),frozenset('ptk'),frozenset('bd'),
			frozenset('ptbdsz'), frozenset('sz'), frozenset('pt'),
			frozenset('pts'), frozenset('ptbd')
		])
		self.assertEqual(classes, self.voice_voiceless.nat_class_set.classes)
		
		
	def test_weight(self):
		self.assertAlmostEqual(self.voice_voiceless.weight('abstakt'), 1.5)
		self.assertAlmostEqual(self.voice_voiceless.weight('skap'), 0.0)
		self.assertAlmostEqual(self.ccvcc.weight('abstakt'), 3)
		self.assertAlmostEqual(self.ccvcc.weight('ppeskutt'), 6)
	
	def test_intersect(self):
		intersect = self.voice_voiceless.intersect(self.ccvcc)
		self.assertAlmostEqual(intersect.weight('abstakt'), 4.5)
	
	def test_error(self):
		arcs = self.voice_voiceless.cat_arcs
		arcs['#'][frozenset(['b'])] = ('#', -1)
		self.assertRaises( ArcLabelError, NaturalClassWFSA.__init__,
			self.voice_voiceless, self.nat_class_set, '#', {'#':-1,'VStop':-1},
			arcs, [1.5] )
		arcs['#'].pop(frozenset(['b']))
		arcs['#'][frozenset('bg')] = ('#', -1)
		self.assertRaises( NoNaturalClassError, NaturalClassWFSA.__init__,
			self.voice_voiceless, self.nat_class_set, '#', {'#':-1,'VStop':-1},
			arcs, [1.5] )
		
	def test_complexity(self):
		pass
	
if __name__ == "__main__":
	unittest.main()
