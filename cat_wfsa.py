from wfsa import *
from param_wfsa import ParametrizedWFSA
from math import log

class CatWFSA(ParametrizedWFSA):
	def __init__(self, alphabet, start, stops, arcs, categories, parameters,
			precision=None, m=None, e=None):
		'''arcs: a dictionary from states to dictionaries
				from letters or frozen sets of categories to 
				state, weight tuples. The categories on an arc 
				label are intersected to get the set of letters
				that may traverse that arc.
		   categories: a dictionary from categories to the set
				of phones belonging to that category.
		'''
		#self.categories[symbol] = set of phonemes
		self.categories = {}
		for cat, phones in categories.items():
			self.categories[cat] = set(phones)
		
		alphabet = set(alphabet)
		
		#convert CatWFSA into an ordinary ParametrizedWFSA
		self.cat_arcs = arcs
		new_arcs = {}
		for state, v in arcs.items():
			new_arcs[state] = {}
			for label, val in v.items():
				if isinstance(label, frozenset):
					phones = alphabet
					for cat in label:
						phones = phones.intersection(self.categories[cat])
					for phone in phones:
						if phone in new_arcs[state].keys():
							#do error checking
							dest0, w0 = new_arcs[state][phone]
							dest1, w1 = val
							if dest0 != dest1 or w0 != w1:
								raise ArcCategoryError( state, phone, 
									dest0, dest1, w0, w1 )
						else:
							new_arcs[state][phone] = val
				else:
					if label in new_arcs[state].keys():
						#do error checking
						dest0, w0 = self.arcs[state][label]
						dest1, w1 = val
						if dest0 != dest1 or w0 != w1:
							raise ArcCategoryError( state, label, 
								dest0, dest1, w0, w1 )
					else:
						new_arcs[state][label] = val
		self.complexity = 1
		ParametrizedWFSA.__init__(self, alphabet, start, stops, new_arcs, parameters)
		self.complexity = self._complex()
	
	def _complex(self):
		states = integer_code_len(len(self.states))
		stops = integer_code_len(len(self.stops))
		param = integer_code_len(len(self.params))
		
		num_arcs = 0
		weight_cost = log(len(self.params)+1, 2)
		labels_len = 0
		for state, out in self.cat_arcs.items():
			for label, dest in out.items():
				num_arcs += 1
				if isinstance(label, frozenset):
					labels_len += integer_code_len(len(label))
					labels_len += len(label)*log(len(self.categories),2)
				else:
					labels_len += log(len(self.alphabet)+1,2)
		state_cost = log(len(self.states),2)
		return num_arcs*(2*state_cost + weight_cost) + labels_len + \
			len(self.stops)*(state_cost + weight_cost) + len(self.params)*self.weight_len\
			+ integer_code_len(num_arcs) + param + states + stops
	
class SingleParamCatWFSA(CatWFSA):
	def __init__( self, alphabet, start, stops, arcs, categories, parameter,
			precision=None, m=None, e=None ):
		CatWFSA.__init__(self, alphabet, start, stops, arcs, categories,
			[parameter], precision, m, e )
		self.complexity -= integer_code_len(1)
	

class ArcCategoryError(Exception):
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