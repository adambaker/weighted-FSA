from math import *
from semiring import *
from arc_set import *
from state import *

'''The classes defined in this module are all deterministic finite state string
 to weight transducers.'''

def integer_code_len( integer, base=2 ):
	integer = abs(integer)
	total = log(2.865064, base)
	current = log(integer, base)
	while( current > 0 ):
		total += current
		current = log(current, base)
	return total

def for_each_arc(arcs, action):
	for source in arcs:
		for label in arcs[source]:
			dest, weight = arcs[source][label]
			action(source, label, dest, weight)

class WeightedFSA(object):
	'''A generic weighted FSA. It is:
			0) A set of symbols, called the alphabet.
			1) A set of states.
			2) A start state and starting weight.
			3) A set of final states and the stopping weight for each final state.
			4) A set of weighted state transition arcs.
			5) A semiring with binary operations for combining arc weights in a path
			and combining path weights.'''
	alt_precision = False
	mantissa = 10
	exp = 5

	def __init__( self, alphabet, start, semiring, start_weight,
				stops=None, arcs=None, arcsets = None, states = None,
				precision = None, m = None,e = None, zero = False):
		'''
		alphabet: the charecter alphabet used for this wfsa. 
		start: the start state. states are intended to be picked out by strings
		semiring: defines a sum, times, zero, and one element. Times is used to
			aggregate weights when traveling a path. Zero, One, and plus are used
			in weight pushing and may be used by subclasses.
		start_weight: the initial weight
		stops: dict(state -> stop weight)
		arcs: dict(state -> dict(letter -> pair(destination state, arc weight)))
		arcsets: a collection of ArcSets
		states: a collection of States. Either states must be specified, or 
			stops and one of arcs or arcsets must be specified.
		precision: boolean. If false, standard 64 bit floats are used as the arc
		weights. Otherwise, arc weights are represented according to m and e.
		'''
		self.start = (start, semiring(start_weight))
		self.semiring = semiring
		self.alphabet = frozenset(alphabet)
		
		if stops is not None:
			if arcs is not None:
				arcsets = self.__arcsets(arcs)
			self._change_precision_and_semiring(arcsets, precision, zero, m, e)
			states = []
			for state_name in set([x.source for x in arcsets]):
				new_state = State(
					state_name, self.alphabet, 
					semiring(stops[state_name]) if state_name in stops else semiring.zero,
					[x for x in arcsets if x.source == state_name]
				)
				states.append(new_state)
			
		self.__states = {}
		for state in states:
			self.__states[state.name] = state
		
		self.state_names = frozenset(self.__states.keys())
		self.trim()
		if not hasattr(self, 'complexity'):
			self.complexity = self._complexity()
						
	def _change_precision_and_semiring(self, arcsets, precision, zero, m, e):
		if precision is None:
			precision = WeightedFSA.alt_precision		
		self.alt_precision = precision
		self.zero = zero
		# round off weights if an alternative to 64 bit floats is used
		if precision:
			if m is None:
				m = WeightedFSA.mantissa
			if e is None:
				e = WeightedFSA.exp	
			self.mantissa = m
			self.exp = e
			self.weight_len = 1 + m + e
			for arcset in arcsets.values():
				arcset.weight = self.round(float(arcset.weight))
		else:
			self.weight_len = 64		
		if zero:
			self.weight_len += 1
		
		for arcset in arcsets:
			arcset.weight = self.semiring(arcset.weight)
		
	def round(self, number):
		m, e = frexp(number)
		if e > 2**(self.exp-1):
			raise WFSAOverflowError(self.mantissa, self.exp, number)
		if e < -(2**(self.exp-1)):
			return 0.0
		m = m*(2**(self.mantissa))
		m = round(m)/float(2**self.mantissa)
		return ldexp(m, e)
	
	def transition(self, state, letter):
		if state not in self.__states:
			return (None, self.semiring.zero)
		ret =  self.__states[state].transition(letter)
		if ret == (None, None):
			return (None, self.semiring.zero)
		else:
			return ret
	
	def stop_weight(self, state_name):
		return self.__states[state_name].stop()
	
	def all_transitions(self):
		for state in self.__states.values():
			for letter in self.alphabet:
				dest, weight = state.transition(letter)
				if dest is not None:
					yield (state.name, letter, dest, weight)
	
	def __remove_other_arcs(self, arcs):
		for state in arcs:
			dest_and_weight = arcs[state].pop('_other', None)
			if dest_and_weight is not None:
				out_letters = set( [x for x in arcs[state].keys()] )
				add_letters = self.alphabet.difference(out_letters)
				for letter in add_letters:
					arcs[state][letter] = dest_and_weight

	def __arcsets( self, arcs ):
		self.__remove_other_arcs(arcs)
		arcsets = {}
		def add_arc_to_arcset(source, label, dest, weight):
			if (source, dest, weight) not in arcsets:
				arcsets[(source,dest,weight)] = []
			arcsets[(source,dest,weight)].append(label)
		
		for_each_arc( arcs, add_arc_to_arcset )
		
		ret = []
		for (source, dest, weight), letters in arcsets.items():
			ret.append(ArcSet(source, dest, letters, weight))
		return ret

	def intersect(self, wfsa ):
		states = []
		
		for state1 in self.__states.values():
			for state2 in wfsa.__states.values():
				states.append( state1.combine(state2) )
		
		start = self.start[0]+'%'+wfsa.start[0]
		start_weight = self.start[1]*wfsa.start[1]
		
		import param_wfsa
		if isinstance(self, MultWFSA):
			new_wfsa = MultWFSA( self.alphabet, start, None, None, states=states)
		elif isinstance(self, param_wfsa.ParametrizedWFSA) \
				and  isinstance(wfsa, param_wfsa.ParametrizedWFSA):
			new_wfsa = param_wfsa.ParametrizedWFSA( self.alphabet, start, states,
				self._params+wfsa._params, precision=False )
		elif isinstance(self, LogWFSA):
			new_wfsa = LogWFSA( self.alphabet, start, None, None, precision=False, 
				states=states )
		else:
			new_wfsa = WeightedFSA( states, start, self.semiring, start_weight,
				precision=False, states=states )
		new_wfsa.complexity = self.complexity + wfsa.complexity
		return new_wfsa
	
	def trim(self):
		reachable = set([self.start[0]])
		
		next = set([])
		for letter in self.alphabet:
			next_state, weight = self.__states[self.start[0]].transition(letter)
			if next_state is not None and weight!= self.semiring.zero:
				next.add(next_state)
		
		while not next.issubset(reachable):
			reachable = reachable.union(next)
			next = set([])
			for state in reachable:
				for letter in self.alphabet:
					next_state, weight = self.__states[state].transition(letter)
					if next_state is not None and weight!= self.semiring.zero:
						next.add(next_state)
		
		end_reachable = set([x.name for x in self.__states.values() 
				if x.stop() != self.semiring.zero])
		current = set(self.__states.keys())
		next = current.difference(end_reachable)
		
		def can_reach_end( state ):
			for letter in self.alphabet:
				dest, weight = state.transition(letter)
				if dest in end_reachable and weight is not None and \
						weight != self.semiring.zero:
					return True
			return False
				
		while next != current and len(next) != 0:
			current = next
			next = set([])
		
			for name in current:
				state = self.__states[name]
				if can_reach_end(state):
					end_reachable.add(name)
				else:
					next.add(name)
			
		#remove unreachable states and associated arcs
		reachable.intersection_update(end_reachable)
		for name in self.state_names.difference(reachable):
			self.__states.pop(name)
			for state in reachable:
				self.__states[state].prune_transitions(name)
		self.state_names = frozenset(self.__states.keys())
	
	def all_pairs_shortest(self):
		'''This is the Gen-All-Pairs algorithm from
		http://www.cs.nyu.edu/~mohri/postscript/hwa.pdf
		For each pair of states in the machine, it finds the sum of all
		the possible paths starting from the first state and ending in 
		the second.
		The results are used for fsa weight pushing.'''

		d = {} #d[s1][s2] is the shortest distance between states s1 and s2
		for s1 in self.__states:
			d[s1] = {}
			for s2 in self.__states:
				d[s1][s2] = self.semiring.zero
				for letter in self.alphabet:
					dest, weight = self.__states[s1].transition(letter)
					if dest == s2:
						d[s1][s2] += weight
		for s0 in self.__states:
			for s1 in self.__states:
				if s1 == s0: continue
				for s2 in self.__states:
					if s2 == s0: continue
					d[s1][s2] += d[s1][s0]*d[s0][s0].star*d[s0][s2]
			for s1 in self.__states:
				if s1 == s0: continue
				d[s0][s1] = d[s0][s0].star*d[s0][s1]
				d[s1][s0] = d[s1][s0]*d[s0][s0].star
			d[s0][s0] = d[s0][s0].star
		return d
	
	def push_weight(self):		
		distance = self.all_pairs_shortest()
		finish = {}
		for state0 in distance:
			sum = self.semiring.zero
			for state1 in distance[state0]:
				sum += distance[state0][state1]*self.stop_weight(state1) 
			finish[state0] = sum
		
		transitions = {}
		def new_weights( source, letter, dest, weight ):
			new_weight = (self.semiring.one/finish[source])*weight*finish[dest]
			if source not in transitions:
				transitions[source] = {}
			transitions[source][letter] = (dest, new_weight)
		
		for source, letter, dest, weight in self.all_transitions():
			new_weights(source, letter, dest, weight)
		
		for state in transitions:
			self.__states[state] = State(
					state, self.alphabet, 
					self.__states[state].stop()/finish[state],
					transitions = transitions[state]
				)
		self.start = (self.start[0], self.start[1]*finish[self.start[0]])
	
	def weight(self, word, bound_strip = True):
		if word[0]=='#' and word[-1]=='#' and bound_strip:
			word = word[1:-1]
		name, weight = self.start
		path = [name]
		for letter in word:
			state = self.__states[name]
			name, w =  state.transition(letter)
			if name is None:
				print word
				print path
				return self.semiring.zero
			weight *= w
			path.append(name)
		weight *= self.__states[name].stop()
		return weight
	
	def print_model(self):
		print self.alphabet
		print self.state_names
		print self.start
		for state in self.__states.values():
			print state.name, state.stop()
			for letter in self.alphabet:
				print state.name, letter, '->', self.transition(state.name, letter)
	
	def _complexity(self):
		num_states = len(self.__states)
		try:
			num_arcs = sum(map(lambda s:s.num_arcs(), self.__states.values()))
		except TypeError:
			return None
		complexity = integer_code_len(num_states) + integer_code_len(num_arcs)
		weight_encoding = \
			lambda x: 1 if float(x) == 0.0 and self.zero else self.weight_len
		for state in self.__states.values():
			complexity += state.complexity(num_states, weight_encoding)
		return complexity
	

class MultWFSA( WeightedFSA ):
	def __init__( self, alphabet, start, stops, arcs, zero=False, states = None ):
		if states is None:
			WeightedFSA.__init__(self, alphabet, start, Probability,
				1.0, stops, arcs, precision=False, zero=zero )
		else:
			WeightedFSA.__init__(self, alphabet, start, Probability, 1.0,
				zero = zero, states = states)
	
	def norm_constant( self, delta=0.000000000001, max_iterations=700 ):
		'''Calculates the total weight for all possible paths through the 
		machine, returning a tuple with the total weight, and True if the
		weight has converged or False if it had not.'''
		
		#This is a dynamic algorithm. Every iteration, it keeps track
		#of all the weight that has exited the machine so far, and all
		#the weight of all the paths that land in each state on that
		#iteration. Carried out to many itterations, this approximates
		#the total weight assigned to the infinite number of paths through
		#the machine, and can be used as a normalizing constant to derive  
		#a probability for a word from the weights the machine assigns to
		#that word.
		
		states = {self.start[0]: Probability.one} 	#maps states to the weight of all
										#path that ends on that state at that
										#iteration
		
		stop = Probability.zero		#total weight of all paths that have stopped
		for i in range(max_iterations):
			next = {}			
			next_stop = Probability.zero	#the ammount to add to stop at 
											#the end of this iteration
			stop_updated = False
			
			for name in states:
				if self.stop_weight(name) != Probability.zero:
					next_stop += states[name]*self.stop_weight(name)
					stop_updated = True
				
			def update_next_state_weights(source, letter, dest, weight):
				if source in states:
					if dest not in next:
						next[dest] = states[source]*weight
					else:
						next[dest] += states[source]*weight
			
			for source, letter, dest, weight in self.all_transitions():
				update_next_state_weights(source, letter, dest, weight)
			
			if stop_updated:
				if float(next_stop) < delta:
					return float(stop + next_stop), True
			states = next
			stop += next_stop
		print max_iterations, 'iterations reached without convergence'
		return float(states[end_state]), False
	
	def log_wfsa(self, base=2):
		states = []
		for state in self.__states.values():
			states.append(state.change_semiring(lambda x:x.to_tropical(base)))
		ret = LogWFSA( self.alphabet, self.start[0], None, None, states=states )
		ret.complexity = self.complexity
		return ret
	

class LogWFSA( WeightedFSA ):
	def __init__( self, alphabet, start, stops, arcs, precision = None, 
			m=None, e=None, zero = False, states=None ):
		if states is None:
			WeightedFSA.__init__(self, alphabet, start, Tropical, 0.0, stops, arcs, 
				precision, m, e, zero)
		else:
			WeightedFSA.__init__(self, alphabet, start, Tropical, 0.0,
				precision=precision, m=m, e=e, zero=zero, states=states)
	
	def mult_wfsa( self, base=2 ):
		states = []
		for state in self.__states.values():
			states.append(state.change_semiring(lambda x:x.to_probability(base)))
		ret = LogWFSA( self.alphabet, self.start[0], None, None, states=states )
		ret.complexity = self.complexity
		return ret
	

class WFSAOverflowError(Exception):
	def __init__( mantissa, exp, number ):
		self.mantissa = mantissa
		self.exp = exp
		self.number = number
	
	def __str__(self):
		return str(self.exp)+' bit exponent is too small to represent '\
			+str(self.number)+'.'
