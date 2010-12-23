import unittest
import fsa.wfsa as wfsa
from math import log, ceil
from fsa.semiring import *


class TestMultFSA( unittest.TestCase ):
	def setUp(self):		
		self.stop = {'0':0.3, '1':0.2}
		self.arcs = {
			'$':{'a':('0',0.3), 'b':('1',0.7)},
			'0':{'a':('0',0.5), 'b':('1', 0.2)},
			'1':{'a':('0',0.7), 'b':('1', 0.1)}
		}
		self.arcs2 = {
			'$':{'a':('0',0.3), 'b':('1',0.7)},
			'0':{'a':('1',0.4), 'b':('0', 0.4)},
			'1':{'a':('0',0.5), 'b':('1', 0.3)}
		}
		self.arcs_other = {
			'$':{'a':('0',0.3), 'b':('1',0.7),'_other':('$',1.0)},
			'0':{'a':('0',0.5), 'b':('1', 0.2),'_other':('$', 1.0)},
			'1':{'a':('0',0.7), 'b':('1', 0.1),'_other':('$', 1.0)}
		}
		self.fsa1 = wfsa.MultWFSA( 'ab', '$', self.stop, self.arcs )
		self.fsa2 = wfsa.MultWFSA( 'ab', '$', self.stop, self.arcs2 )
		self.other = wfsa.MultWFSA( 'pernicious ab', '$', self.stop, self.arcs_other )
	
	def test_other(self):
		ab_weight = float(self.fsa1.weight('ab'))
		self.assertAlmostEqual(0.012, ab_weight)
		cab_weight = float(self.other.weight('cab'))
		self.assertAlmostEqual(ab_weight, cab_weight)
		bca_weight = float(self.other.weight('bca'))
		self.assertAlmostEqual(0.063, bca_weight)
		self.assertAlmostEqual(ab_weight, float(self.other.weight('pernicious ab')))
		
	def test_intersect(self):
		stops = {'0%0':Probability(0.09),'0%1':Probability(0.06),
				 '1%0':Probability(0.06),'1%1':Probability(0.04)} 
		arcs = \
		{
			'$%$':{'a':('0%0', Probability(0.09)), 'b':('1%1', Probability(0.49))},
			'0%0':{'a':('0%1', Probability(0.2)), 'b':('1%0', Probability(0.08))},
			'0%1':{'a':('0%0', Probability(0.25)), 'b':('1%1', Probability(0.06))},
			'1%0':{'a':('0%1', Probability(0.28)), 'b':('1%0', Probability(0.04))},
			'1%1':{'a':('0%0', Probability(0.35)), 'b':('1%1', Probability(0.03))}
		}
		intersect = self.fsa1.intersect(self.fsa2)
		self.assertEquals(intersect.start, ('$%$', Probability(1.0)))
		
		def test_arc( source, letter, dest, weight ):
			dest0, weight0 = intersect.transition(source, letter)
			self.assertEquals(dest, dest0)
			self.assertAlmostEquals(float(weight), float(weight0))
		
		wfsa.for_each_arc( arcs, test_arc )
		for state in stops:
			self.assertAlmostEquals( float(intersect.stop_weight(state)), 
					float(stops[state]) )
	
	def test_norm(self):
		self.assertAlmostEqual(self.fsa1.norm_constant()[0], 1.0)
		new = (self.arcs['0']['b'][0], self.arcs['0']['b'][1]+0.1)
		self.arcs['0']['b'] = new
		self.fsa1 = wfsa.MultWFSA('ab', '$', self.stop, self.arcs)
		self.assert_(self.fsa1.norm_constant()[0] > 1.0)
		new = (self.arcs['0']['b'][0], self.arcs['0']['b'][1]-0.2)
		self.arcs['0']['b'] = new
		self.fsa1 = wfsa.MultWFSA('ab', '$', self.stop, self.arcs)
		self.assert_(self.fsa1.norm_constant()[0] < 1.0)
	
	def test_trim(self):
		self.arcs['0'].pop('b')
		self.arcs['$'].pop('b')
		
		self.assertEqual(
			wfsa.MultWFSA('ab', '$', self.stop, self.arcs).state_names,
			set(['0','$'])
		)
		
		self.arcs2['1']['a'] = ('1', 0.5)
		self.stop.pop('1')
		
		self.assertEqual(
			wfsa.MultWFSA('ab', '$', self.stop, self.arcs).state_names,
			set(['0','$'])
		)
	
	def test_complexity(self):
		#upper bound the complexity value: encode number of states, number of arcs
		#each arc, each stop (with no stop for a state as 0), number of stops
		state = wfsa.integer_code_len(3)
		arc = wfsa.integer_code_len(6)
		stop = wfsa.integer_code_len(2)
		single_arc = 3*log(3,2)+64
		single_stop = log(3,2)+64
		complexity = 6*single_arc + 3*single_stop + arc + state + stop
		self.assert_( self.fsa1.complexity <= complexity )
				
		self.arcs['0']['a'] = ('0', 0.0)
		fsa3 = wfsa.MultWFSA( self.fsa1.alphabet, self.fsa1.start[0], self.stop, 
			self.arcs, zero = True )
		self.assert_( fsa3.complexity, self.fsa1.complexity )
	
	def test_int_code(self):
		self.assertAlmostEquals(wfsa.integer_code_len(3), 3.767978574539)
		self.assertAlmostEquals(wfsa.integer_code_len(4), 4.518567366364)
	
	def test_all_pairs(self):
		self.arcs['$']['a']=('0',Probability(0.5))
		self.arcs['1']['b'] = ('1',Probability(0.5))
		self.fsa1 = wfsa.MultWFSA('ab', '$', self.stop, self.arcs )
		d = self.fsa1.all_pairs_shortest()
		
		finish = Probability.zero#calculating the total distance from the start state
					#to all final states, times the weight of exiting at 
					#each final state. This value should be the same as
					#the FSA's norm constant.
		start = self.fsa1.start[0]
		finish = self.fsa1.stop_weight(start)
		for state_name in d[start]:
			finish += d[start][state_name]*self.fsa1.stop_weight(state_name)
		norm = self.fsa1.norm_constant()
		self.assert_(norm[1])
		self.assertAlmostEqual(float(finish), float(norm[0]), 5)
	
	def test_weight_push(self):
		self.arcs['$']['a']=('0',Probability(0.5))
		self.arcs['1']['b'] = ('1',Probability(0.5))
		self.fsa1 = wfsa.MultWFSA('ab', '$', self.stop, self.arcs )
		ab_weight = self.fsa1.weight('ab')
		abaaab_weight = self.fsa1.weight('abaaab')
		norm = self.fsa1.norm_constant()[0]
		self.fsa1.push_weight()
		self.assertAlmostEqual(float(self.fsa1.start[1]), norm, 5)
		for state in self.fsa1.state_names:
			sum = self.fsa1.stop_weight(state)
			for letter in self.fsa1.alphabet:
				dest, weight = self.fsa1.transition(state, letter)
				sum += weight
			self.assertAlmostEqual( float(sum), 1.0 )
		self.assertAlmostEqual( float(self.fsa1.weight('ab')), float(ab_weight) )
		self.assertAlmostEqual( float(self.fsa1.weight('abaaab')), float(abaaab_weight) )

if __name__ == "__main__":
	unittest.main()
