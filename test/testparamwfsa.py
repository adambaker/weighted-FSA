import unittest
import param_wfsa as fsa
from arc_set import ArcSet
from semiring import *

class TestParamFSA( unittest.TestCase ):
	def setUp(self):
		self.parameter1 = [2, 3, 10]
		self.parameter2 = [1, 6, 4]
		self.arcsetStart1 = ArcSet('$', '1', 'ab', -1)
		self.arcsetStart2 = ArcSet('$', '2', 'cd', 0)
		self.arcset11 = ArcSet('1', '1', 'bc', 2)
		self.arcset12 = ArcSet('1', '2', 'ad', -1)
		self.arcset21 = ArcSet('2', '1', 'abc', 1)
		self.arcset22 = ArcSet('2', '1', 'd', -1)
		self.stops = {'1':-1, '2':0}
		
		self.stateset1 = fsa.parametrized_states('abcd', 
				[
					self.arcsetStart1, self.arcsetStart2, self.arcset11, 
					self.arcset12, self.arcset21, self.arcset22
				],
				self.stops, self.parameter1)
		self.fsa1 = fsa.ParametrizedWFSA('abcd', '$', self.stateset1, self.parameter1)
		self.stateset2 = fsa.parametrized_states('abcd', 
				[
					self.arcsetStart1, self.arcsetStart2, self.arcset11, 
					self.arcset12, self.arcset21
				],
				self.stops, self.parameter2)
		self.fsa2 = fsa.ParametrizedWFSA('abcd', '$', self.stateset2, self.parameter2)
	
	def test_weight(self):
		self.assertAlmostEquals( float(self.fsa1.weight('acd')), 12 )
		self.assertAlmostEquals( float(self.fsa1.weight('ccda')), 8 )
		self.fsa1.set_parameter( 0, 4 )
		self.assertAlmostEquals( float(self.fsa1.weight('acd')), 14 )
		self.assertAlmostEquals( float(self.fsa1.weight('ccda')), 10 )
		self.fsa1.set_parameter( 1, 5 )
		self.assertAlmostEquals( float(self.fsa1.weight('acd')), 14 )
		self.assertAlmostEquals( float(self.fsa1.weight('ccda')), 14 )
		self.fsa1.set_parameter( 2, 7 )
		self.assertAlmostEquals( float(self.fsa1.weight('acd')), 11 )
		self.assertAlmostEquals( float(self.fsa1.weight('ccda')), 14 )
	
	def test_intersect(self):
		fsa3 = self.fsa1.intersect(self.fsa2)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 29 )
		fsa3.set_parameter(0, 1)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 27 )
		fsa3.set_parameter(3, 3)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 31 )
		fsa3.set_parameter(1, 1)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 29 )
		fsa3.set_parameter(4, 3)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 26 )
		fsa3.set_parameter(2, 5)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 21 )
		fsa3.set_parameter(5, 5)
		self.assertAlmostEquals( float(fsa3.weight('dcba')), 22 )
	

if __name__ == "__main__":
	unittest.main()