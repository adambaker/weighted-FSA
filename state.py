from arc_set import *
from math import log
import semiring

class State(object):
	'''
	>>> set1 = ArcSet('1', '2', 'aeiou', 4)
	>>> set2 = ArcSet('1', '3', 'fghj', 12)
	>>> set3 = ArcSet('1', '5', 'bcd', 3)
	>>> state = State('1', 'abcdefghijou', 2, [set1, set2, set3])
	>>> state.transition('a')
	('2', 4)
	>>> state.transition('b')
	('5', 3)
	>>> state.transition('c')
	('5', 3)
	>>> state.transition('f')
	('3', 12)
	>>> state2 = State('1', 'bcdefghijou', 3, [set1, set2])
	>>> state2.transition('a')
	('2', 4)
	>>> state2.transition('b')
	(None, None)
	>>> state3 = state.combine(state2)
	>>> state3.transition('a')
	('2%2', 16)
	>>> state3.transition('b')
	(None, None)
	>>> state3.stop()
	6
	>>> state4 = state.change_semiring(lambda x:semiring.Tropical(x))
	>>> state4.transition('a')
	('2', 4 in a tropical semiring.)
	>>> state4.stop()
	2 in a tropical semiring.
	>>> state5 = state4.change_semiring(lambda x:x.to_probability())
	>>> state5.transition('a')
	('2', 0.0625 in a probability semiring.)
	>>> state5.stop()
	0.25 in a probability semiring.
	'''
	def __init__(self, name, alphabet, stop_weight, arcsets=None, transitions=None):
		self.name = name
		self._alphabet = frozenset(alphabet)
		self._stop_weight = stop_weight
		if transitions is not None:
			self._transitions = transitions
			self._arcs = None
		else:
			if arcsets is None:
				raise TypeError('Must specify either arcsets or transitions')
			self.__arcsets = arcsets
			self._transitions = {}
			for arcset in arcsets:
				if arcset.source != name:
					raise IllegalState(
						"The State's name must be identical to the source " +
						"for each ArcSet. " + name + ' != ' + arcset.source
					)
				for letter in arcset:
					self._transitions[letter] = (arcset.dest, arcset.weight)
			self._arcs = self._optimal_arcs(arcsets)
	
	def stop(self):
		return self._stop_weight
		
	def transition(self, letter):
		return self._transitions[letter] \
			if letter in self._transitions else (None, None)
		
	def combine( self, other ):
		name = self.name + '%' + other.name
		stop_weight = self.stop()*other.stop()
		transitions = {}
		for letter, (dest1, w1) in self._transitions.items():
			dest2, w2 = other.transition(letter)
			if dest2 is not None:
				transitions[letter] = (dest1+'%'+dest2, w1*w2)
		return State(name, self._alphabet, stop_weight, transitions=transitions)
		
	def _optimal_arcs(self, arcsets):
		largest = set([])
		letters_covered = frozenset([])
		for arcset in arcsets:
			letters_covered = letters_covered.union(arcset.letters)
			if len(largest) < len(arcset):
				largest = arcset
		arcs = {}
		for arcset in arcsets:
			if arcset == largest and letters_covered == self._alphabet:
				arcs['_other'] = (arcset.dest, arcset.weight)
			else:
				for letter in arcset:
					arcs[letter] = (arcset.dest, arcset.weight)
		self._validate_arcsets(arcsets)
		return arcs
	
	def _validate_arcsets(self, arcsets):
		letters = set([])
		for arcset in arcsets:
			if len(arcset.letters.intersection(letters)) != 0:
				raise IllegalState('The letter '+letter+
									' appears in more than one ArcSet.' 
				)
			letters.update(arcset.letters)
	
	def num_arcs(self):
		if self._arcs == None:
			return None
		else:
			return len(self._arcs)
	
	def complexity(self, num_states, weight_encoding, base=2):
		if self._arcs == None:
			return None
		num_labels = len(self._alphabet)+1
		state_enc = log(num_states, base)
		labels_encoded = 0
		tot_complexity = weight_encoding(self._stop_weight)
		for letter, (dest, weight) in self._arcs.items():
			tot_complexity += 2*state_enc+weight_encoding(weight)
			tot_complexity += log(num_labels-labels_encoded, base)
			labels_encoded += 1
		return tot_complexity

	def prune_transitions(self, dest):
		arc_removed = False
		for label in self._transitions:
			if self._transitions[label][0] == dest:
				self._transitions.pop(label)
				arc_removed = True
				if self._arcs is not None:
					self._arcs.pop(label)
		if arc_removed and self._arcs is not None and '_other' in self._arcs:
			self._arcs.pop('_other')
			for letter in self._alphabet.difference(set(self._arcs.keys())):
				self._arcs[letter] = self.__transitions[letter]
	
	def normalize(self):
		total = sum([x[1] for x in self.__transitions.values()])
		total += self._stop_weight

		trans = {}
		for letter, (dest, weight) in self.__transitions.items():
			trans[letter] = (dest, weight/total)
		return State(self.name, self._alphabet, 
				self._stop_weight/total, transitions=trans)
	
	def change_semiring(self, change):
		transitions = {}
		for letter, (dest, weight) in self._transitions.items():
			transitions[letter] = (dest, change(weight))
		return State(self.name, self._alphabet, change(self._stop_weight),
			transitions=transitions)
	

class ParametrizedState(State):
	'''
	A State whose arc weights can be tied and edited dynamically.
	
	>>> from semiring import Tropical
	>>> alphabet = set('abcdef')
	>>> trans1 = {'a': ('1', Tropical(3.0)), 'b':('1', Tropical(2.0)),
	... 'c':('2',Tropical(3.0)), 'd':('3',Tropical(2.0)), 
	... 'e':('1',Tropical(3.0)), 'f':('2',Tropical(1.0)) }
	>>> trans_map1 = [set('f'), set('db'), set('ace')]
	>>> trans2 = {'a': ('2', Tropical(4.0)), 'b':('1', Tropical(3.0)),
	... 'c':('1', Tropical(4.0)), 'd':('1', Tropical(4.0)),
	... 'e':('2', Tropical(4.0)), 'f':('2', Tropical(4.0))}
	>>> trans_map2 = [set(['_stop','b']), set('acdef')]
	>>> state1 = ParametrizedState('1', alphabet, Tropical.one, trans_map1,
	... transitions=trans1 )
	>>> state2 = ParametrizedState('2', alphabet, Tropical(3.0), trans_map2,
	... transitions=trans2 )
	>>> state1.transition('a')
	('1', 3.0 in a tropical semiring.)
	>>> state1.change_parameter(2, Tropical(6.0), Tropical(3.0))
	>>> state1.transition('a')
	('1', 6.0 in a tropical semiring.)
	>>> state2.transition('b')
	('1', 3.0 in a tropical semiring.)
	>>> state2.stop()
	3.0 in a tropical semiring.
	>>> state2.change_parameter(0, Tropical(5.0), Tropical(3.0))
	>>> state2.transition('b')
	('1', 5.0 in a tropical semiring.)
	>>> state2.stop()
	5.0 in a tropical semiring.
	>>> state3 = state1.combine(state2)
	>>> state3.transition('f')
	('2%2', 5.0 in a tropical semiring.)
	>>> state3.change_parameter(4, Tropical(8.0), Tropical(4.0))
	>>> state3.transition('f')[1]
	9.0 in a tropical semiring.
	>>> state3.transition('a')[1]
	14.0 in a tropical semiring.
	>>> state3.change_parameter(1, Tropical(0.0), Tropical(2.0))
	>>> state3.transition('d')[1]
	8.0 in a tropical semiring.
	>>> state3.transition('b')[1]
	5.0 in a tropical semiring.
	'''
	def __init__(self, name, alphabet, stop_weight, parameter_map,
			arcsets=None, transitions=None):
		'''
		All parameters but parameter_map function exactly as in State.
		parameter_map is a list of collections of transitions whose weights are tied
			to the parameter at that collection's index. The special transition label
			"_stop" indicates the stopping weight for this state is tied to the 
			indexed parameter.
		'''
		State.__init__(self, name, alphabet, stop_weight, arcsets, transitions)
		self._parameter_map = parameter_map
		
	def uses_parameter(self, index):
		return len(self._parameter_map[index]) > 0
		
	def change_parameter(self, index, value, old_value):
		for transition in self._parameter_map[index]:
			if transition == '_stop':
				self._stop_weight /= old_value
				self._stop_weight *= value
			else:
				try:
					self._transitions[transition] = \
						(
							self._transitions[transition][0], 
							self._transitions[transition][1]*value/old_value
						)
				except:
					print self.name
	
	def combine( self, other ):
		state = State.combine(self, other)
		new_param_map = []
		new_param_map.extend(self._parameter_map)
		new_param_map.extend(other._parameter_map)
		return ParametrizedState(state.name, self._alphabet, state.stop(),
				new_param_map, transitions=state._transitions)
	
	def prune_transitions(self, dest):
		State.prune_transitions(self, dest)
		for param_point in self._parameter_map:
			for label in frozenset(param_point):
				if label != '_stop' and self.transition(label)[0] is None:
					param_point.remove(label)
	
class NaturalClassState(ParametrizedState):
	'''
	>>> from nat_class_set import NaturalClassSet
	>>> from semiring import Tropical
	>>> coronal = frozenset('tpbdrzlnm')
	>>> stop = frozenset('bd')
	>>> nasal = frozenset('nm')
	>>> classes = NaturalClassSet('ptkbdgrszlmn', set([coronal, nasal, stop]))
	>>> arcset = NatClassArcSet('1', '2', 'bdrzlmn', Tropical(3.0) , classes)
	>>> arcset2 = NatClassArcSet('1', '1', 'ptkgs', Tropical.one, classes)
	>>> param_map = [set('bdrzlmn')]
	>>> state = NaturalClassState('1', classes, Tropical.one, param_map,
	... [arcset, arcset2])
	>>> state.transition('b')
	('2', 3.0 in a tropical semiring.)
	>>> state.transition('p')
	('1', 0.0 in a tropical semiring.)
	>>> len(state._arcs)
	6
	>>> state._arcs['_other']
	('1', 0.0 in a tropical semiring.)
	>>> state._arcs[stop]
	('2', 3.0 in a tropical semiring.)
	>>> state._arcs[nasal]
	('2', 3.0 in a tropical semiring.)
	>>> state._arcs[frozenset('z')]
	('2', 3.0 in a tropical semiring.)
	>>> state._arcs[frozenset('l')]
	('2', 3.0 in a tropical semiring.)
	>>> state._arcs[frozenset('r')]
	('2', 3.0 in a tropical semiring.)
	'''
	def __init__(self, name, nat_class_set, stop_weight, parameter_map, 
					arcsets=None, transitions=None):
		self.nat_class_set = nat_class_set
		ParametrizedState.__init__(self, name, nat_class_set.alphabet, 
									stop_weight, parameter_map, arcsets,
									transitions)
	
	def complexity(self, num_states, weight_encoding, base=2):
		if self._arcs == None:
			return None
		state_enc = log(num_states, base)
		param_enc = log(len(self._parameter_map)+1, base)
		labels_size_left = self.nat_class_set.labels_len
		tot_complexity = param_enc
		for label, (dest, weight) in self._arcs.items():
			tot_complexity += 2*state_enc+param_enc
			tot_complexity += log(len(label)/labes_size_left, base)
			labels_size_left -= len(label)
		return tot_complexity

	
	def _optimal_arcs(self, arcsets):
		norm = float(self.nat_class_set.labels_len)
		lowest_prob = 1.0
		lowest_arcset = None
		letters_covered = frozenset([])
		for arcset in arcsets:
			letters_covered = letters_covered.union(arcset.letters)
			curr_prob = 1.0
			for label in arcset.labels:
				curr_prob *= len(label)/norm
			if curr_prob < lowest_prob:
				lowest_prob = curr_prob
				lowest_arcset = arcset
		arcs = {}
		for arcset in arcsets:
			if arcset is lowest_arcset and letters_covered == self._alphabet:
				arcs['_other'] = (arcset.dest, arcset.weight)
			else:
				for label in arcset.labels:
					arcs[label] = (arcset.dest, arcset.weight)
		self._validate_arcsets(arcsets)
		return arcs
			

class IllegalState(Exception):
	def __init__(self, message):
		self.message = message
	
	def __str__(self):
		return self.message
	

if __name__ == '__main__':
	import doctest
	doctest.testmod()
