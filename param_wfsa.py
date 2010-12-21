from wfsa import *
from semiring import *
from arc_set import ArcSet

def parametrized_states( alphabet, arcsets, stops, parameters, default=Tropical.one ):
	state_names = set([x.source for x in arcsets])
	
	#constructing a parameter map for each state and changing the arcset
	#weight to the weight on the parameter it points to
	parameter_maps = {}
	converted_arcsets = []

	for a in arcsets:
		#make a copy so that the arcset can be reused in specifying another machine
		#if appropriate.
		arcset = ArcSet(a.source, a.dest, a.letters, a.weight)
		if arcset.source not in parameter_maps:
			parameter_maps[arcset.source] = []
			for i in range(len(parameters)):
				parameter_maps[arcset.source].append(set([]))
		if arcset.weight >= 0:
			parameter_maps[arcset.source][arcset.weight].update(arcset.letters)
			arcset.weight = Tropical(parameters[arcset.weight])
		else:
			arcset.weight = default
		converted_arcsets.append(arcset)
	
	states = []
	converted_stops = {}
	for name in state_names:
		if name in stops:
			if stops[name] >= 0:
				parameter_maps[name][stops[name]].add('_stop')
				converted_stops[name] = Tropical(parameters[stops[name]])
			else:
				converted_stops[name] = default
		else:
			converted_stops[name] = Tropical.zero
		
	for name in state_names:
		states.append( ParametrizedState(name, alphabet, converted_stops[name], 
			parameter_maps[name],[x for x in converted_arcsets if x.source == name]
		))
	return states

class ParametrizedWFSA( LogWFSA ):
	def __init__(self, alphabet, start, states, parameters, default=Tropical.one,
			precision = None, m=None, e=None):
		''' 
			parameters: a list of arc weights
			arcs: dict{source state -> dict{letter->(destination state, parameter number)}}
				the parameter number is the index in the parameter list for the relevant
				weight to assign the transition, or a negative number to use the default
				weight.
		'''
		self._params = parameters
		self.__states_with_parameter = []
		for index in range(len(parameters)):
			self.__states_with_parameter.append(
				[s for s in states if s.uses_parameter(index)]
			)
		self.complexity = 0
		LogWFSA.__init__( self, alphabet, start, None, None, precision, 
				m, e, False, states )
		#self.complexity = self._complexity()
	
	def set_parameter(self, index, value):
		old_value = self._params[index]
		for state in self.__states_with_parameter[index]:
			state.change_parameter(index, Tropical(value), Tropical(old_value))
		self._params[index] = value
		
	def trim(self):
		WeightedFSA.trim(self)
		for index, states in enumerate(self.__states_with_parameter):
			self.__states_with_parameter[index] = \
					filter( lambda s: s.name in self.state_names, states )
					
	
	def _complexity(self):
		states = integer_code_len(len(self.states))
		param = integer_code_len(len(self.params))
		
		num_arcs = 0
		weight_cost = log(len(self.params)+1, 2)
		for state, out in self.arcs.items():
			for letter, dest in out.items():
				num_arcs += 1
		state_cost = log(len(self.states),2)
		return num_arcs*(2*state_cost + log(len(self.alphabet)+1,2) + weight_cost) + \
			len(self.stops)*(state_cost + weight_cost) + len(self.params)*self.weight_len\
			+ integer_code_len(num_arcs) + param + states + stops

