from wfsa import *
from param_wfsa import ParametrizedWFSA
from math import log
	

class NaturalClassWFSA( ParametrizedWFSA ):
	'''
	A ParametrizedWFSA whose arc labels are sets of phones corresponding to
		a natural class.
	classes: A superset of the natural classes used in this WFSA. Each
			natural class is a set of phonemes in that class. The full set of
			natural classes is gotten by closing this set under intersection 
			and adding all the phonemes in alphabet.
	'''
		
	def __init__( self, nat_class_set, start, states, parameters,
			 precision=None, m=None, e=None ):
		'''See: param_wfsa.ParametrizedWFSA
		arcs: dict(start state -> dict(natural class -> 
				pair(destination state, parameter number) ) )
		nat_class_set: A NaturalClassSet, which contains all the natural classes
				in the language, plus the alphabet, plus information useful for
				determining machine complexity.
		'''
		self.nat_class_set = nat_class_set
		self.alphabet = nat_class_set.alphabet
		self.complexity = 0
		ParametrizedWFSA.__init__( self, self.alphabet, start, states, 
				parameters, 0, precision, m, e )
		self.complexity = self._complex()
	
	def _complex(self):
		states = integer_code_len(len(self.states))
		stops = integer_code_len(len(self.stops))
		param = integer_code_len(len(self.params))
		
		weight_cost = log(len(self.params)+1, 2)
		
		#put num_arcs and labels_len in singleton lists so they can be \
		#updated in a callback function.
		num_arcs = [0] 
		labels_len = [0]
		
		def count_arcs_and_sum_codelength_of_labels(source, label, dest, weight):
			num_arcs[0] += 1
			if label == '_other':
				labels_len[0] += -log(
					len(self.alphabet)/self.nat_class_set.labels_len, 
					2
				)
			else:
				labels_len[0] += -log(
					len(label)/self.nat_class_set.labels_len, 
					2
				)
				
		for_each_arc( self.cat_arcs, count_arcs_and_sum_codelength_of_labels )
			
		#unpackage num_arcs and labels_len
		num_arcs = num_arcs[0]
		labels_len = labels_len[0]
		
		state_cost = log(len(self.states),2)
		return num_arcs*(2*state_cost + weight_cost) + labels_len + \
			len(self.stops)*(state_cost + weight_cost) + len(self.params)*self.weight_len\
			+ integer_code_len(num_arcs) + param + states + stops
	
class CountingWFSA( NaturalClassWFSA ):
	'''A NaturalClassWFSA with a single parameter fixed to 1.0. It can be
		used to count how many times a configuration occurs in a string.'''
	def __init__(self, alphabet, start, stops, arcs, classes):
		NaturalClassWFSA.__init__(self, alphabet, start, stops, arcs,
				[1.0], classes, False)
		self.complexity -= integer_code_len(1) #don't need to store # of params
		self.complexity -= self.weight_len #don't need to encode a parameter
	

class NoNaturalClassError(Exception):
	def __init__(self, bad_class):
		self.bad_class = bad_class
	
	def __str__(self):
		return 'The class '+str(self.bad_class)+' is not in the NaturalClassSet.'

class ArcLabelError(Exception):
	def __init__(self, state, letter, dest1, dest2, w1, w2):
		self.state = state
		self.letter = letter
		self.dest1 = dest1
		self.dest2 = dest2
		self.w1 = w1
		self.w2 = w2
	
	def __str__(self):
		return 'Parsing '+self.letter+' from state '\
			+self.state+'leads to both '+self.dest1+\
			' with weight '+self.w1+' and '+self.dest2\
			+' with weight '+self.w2+'.'