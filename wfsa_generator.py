from wfsa import *
from nat_class_wfsa import *

def bigram_plog_wfsa( cond_plogs, zero = False ):
	states = set(['#'])
	start = '#'
	stops = {}
	arcs = {}
	alphabet = set([])
	for bigram in cond_plogs.keys():
		alphabet.add(bigram[0])
		if bigram[1] == '#':
			stops[bigram[0]] = cond_plogs[bigram]
			states.add(bigram[0])
		else:
			if bigram[0] not in arcs.keys():
				arcs[bigram[0]] = {}
			arcs[bigram[0]][bigram[1]] = (bigram[1], cond_plogs[bigram])
			states.add(bigram[0])
			states.add(bigram[1])
	return LogWFSA( alphabet, start, stops, arcs, zero=zero )

def bigram_prob_wfsa( conditionals, zero = False ):
	states = set(['#'])
	start = '#'
	stops = {}
	arcs = {}
	alphabet = set([])
	for bigram in conditionals.keys():
		alphabet.add(bigram[0])
		if bigram[1] == '#':
			stops[bigram[0]] = conditionals[bigram]
			states.add(bigram[0])
		else:
			if bigram[0] not in arcs.keys():
				arcs[bigram[0]] = {}
			arcs[bigram[0]][bigram[1]] = (bigram[1], conditionals[bigram])
			states.add(bigram[0])
			states.add(bigram[1])
			
	return MultWFSA( alphabet, start, stops, arcs, zero=zero )

def v_mi_wfsa( alphabet, v_mi, zero = False ):
	vowels = set([])
	states = set(['#']) 
	start = '#'
	stops = {}
	arcs = {}
	for bigram in v_mi.keys():
		vowels.add(bigram[0])
		#state of having just seen vowel, no following consonant
		states.add(bigram[0])
		#state of having seen vowel and one or more consonants
		states.add(bigram[0]+'CONS')
		if bigram[1] !='#':
			states.add(bigram[1])
			states.add(bigram[1]+'CONS')
			#add V0CONS -> V1 state transitions
			if bigram[0]+'CONS' not in arcs.keys():
				arcs[bigram[0]+'CONS'] = {}
			arcs[bigram[0]+'CONS'][bigram[1]] = (bigram[1], -v_mi[bigram])
		else:
			stops[bigram[0]+'CONS'] = -v_mi[bigram]
	
	alphabet = set(alphabet)
	consonants = alphabet.difference(vowels)
	consonants.discard('#')
		
	if '#' not in arcs.keys():
		arcs['#'] = {}
		
	for consonant in consonants:
		for state in states:
			if state not in arcs.keys():
				arcs[state] = {}
			
			#add V -> VCONS and # -> V transitions
			if state in vowels:
				arcs['#'][state] = (state, 0.0)
				arcs[state][consonant] = (state+'CONS', 0.0)
				
			#add VCONS -> VCONS transition
			else:
				arcs[state][consonant] = (state, 0.0)
				
	#add V -> V transitions
	for first in vowels:
		for second in vowels:
			arcs[first][second] = (second, 0.0)
			
	#add trivial stops
	for state in states:
		if state not in stops.keys():
			stops[state] = 0.0
	
	return LogWFSA( alphabet, start, stops, arcs, zero=zero )

def trigram_stress_wfsa( alphabet, stress, conditionals, zero = False ):
	'''	
	stress: a dictionary from stressed vowels to the vowels' stress values
	
	conditionals: a dictionary from two or three stress values to a 
		positive log conditional probability of the final stress 
		value given the previous value(s). If the key is only two values,
		the first should be the word boundary, '#'

	'''
	stress['#'] = '#'
	vowels = []
	for v in stress.keys():
		vowels.append(v)
	vowels = set(vowels)
	vowels.discard('#')
	alphabet = set(alphabet)
	consonants = alphabet.difference(vowels)
	consonants.discard('#')
	
	states = set(['#'])
	start = '#'
	stops = {}
	arcs = {}
	
	for gram in conditionals.keys():
		#gram is 2- or 3-tuple of stress values
		if len(gram) == 3:
			states.add( gram[0] + ',' + gram[1] )
	
	for state in states:
		arcs[state] = {}
		for c in consonants:
			arcs[state][c] = (state, 0.0)
		
		for v in vowels:
			s = stress[v]
			if state == '#':
				w = conditionals[('#',s)]
				dest = '#,'+s
			else:
				a, b = state.split(',')
				w = conditionals[(a,b,s)]
				dest = b + ',' + s
			arcs[state][v] = (dest, w)
		
		if state == '#':
			stops[state] = conditionals[('#','#')]
		else:
			a, b = state.split(',')
			stops[state] = conditionals[(a,b,'#')]
	
	return LogWFSA( alphabet, start, stops, arcs, zero=zero )

def bigram_stress_wfsa( alphabet, stress, conditionals, zero = False ):
	'''	
	stress: a dictionary from stressed vowels to the vowels' stress values
	
	conditionals: a dictionary from two stress values to a 
		positive log conditional probability of the final stress 
		value given the previous value. 
	'''
	vowels = []
	for v in stress.keys():
		vowels.append(v)
	vowels = set(vowels)
	vowels.discard('#')
	alphabet = set(alphabet)
	alphabet.discard('#')
	consonants = alphabet.difference(vowels)
	
	states = set(['#'])
	start = '#'
	stops = {}
	arcs = {}
	
	for bigram in conditionals.keys():
		states.add( bigram[0] )
	
	for state in states:
		arcs[state] = {}
		for c in consonants:
			arcs[state][c] = (state, 0.0)
		
		for v in vowels:
			s = stress[v]
			w = conditionals[(state,s)]
			arcs[state][v] = (s, w)
		
		stops[state] = conditionals[(state,'#')]
	
	return LogWFSA( alphabet, start, stops, arcs, zero = zero )

def stress_v_mi_wfsa( stress_v_mi, alphabet, zero = False ):
	'''
	stress_v_mi: a dictionary from stressed vowels to the mutual 
		information between the value and it's stress type.
	'''
	states = set(['#'])
	start = '#'
	stops = {'#':0.0}
	arcs = {'#':{}}
	alphabet = set(alphabet)
	alphabet.discard('#')
	
	for letter in alphabet:
		w = 0.0
		if letter in stress_v_mi.keys():
			w = -stress_v_mi[letter]
		arcs['#'][letter] = ('#', w)
	
	return LogWFSA( alphabet, start, stops, arcs, zero = zero )
	
def class_bigram( nat_class_set, class1, class2, param = 1.0 ):
	'''
	>>> consonant = frozenset('ptbdszmn')
	>>> obst = frozenset('ptbdsz')
	>>> voice = frozenset('bdzmn')
	>>> vowel = frozenset('aeiou')
	>>> classes = NaturalClassSet('ptbdszmnaeiou', set([consonant, obst, voice, vowel]))
	>>> obst_voice = class_bigram( classes, obst, voice )
	>>> obst_voice.weight('patzdom')
	2.0
	>>> obst_voice.weight('pazd')
	1.0
	>>> obst_voice.weight('staduot')
	0.0
	'''
	classes = nat_class_set.classes
	class1 = frozenset(class1)
	class2 = frozenset(class2)
	#state mnemonics: '#': no part of pattern processed
	#				  '1': character from class1 processed;
	#					   looking for character from class2
	start = '#'
	stops = {'#':-1, '1':-1}
	arcs = {}
	arcs['#'] = {}
	arcs['#'][class1] = ('1', -1)
	arcs['#']['_other'] = ('#', -1)
	arcs['1'] = {}
	intersect = frozenset(class1.intersection(class2))
	if len(intersect) == 0:
		arcs['1'][class1] = ('1', -1)
		arcs['1'][class2] = ('#', 0)
	else:
		arcs['1'][intersect] = ('1', 0)
		cover_labels = get_covering_labels(class1.difference(intersect), classes)
		for label in cover_labels:
			arcs['1'][label] = ('1', -1)
		cover_labels = get_covering_labels(class2.difference(intersect), classes)
		for label in cover_labels:
			arcs['1'][label] = ('#', 0)
	arcs['1']['_other'] = ('#', -1)
	return NaturalClassWFSA( nat_class_set, start, stops, arcs, [param] )

def class_tier_bigram( nat_class_set, tier, class1, class2, param = 1.0 ):
	'''
	>>> vowel = frozenset('ieaouy')
	>>> consonant = frozenset('stpbdzl')
	>>> front = frozenset('iea')
	>>> back = frozenset('aou')
	>>> alphabet = vowel.union(consonant)
	>>> classes = NaturalClassSet(alphabet, set([vowel, consonant, front, back]))
	>>> front_back = class_tier_bigram( classes, vowel, front, back )
	>>> front_back.weight('stopao')
	1.0
	>>> front_back.weight('itzabsu')
	2.0
	>>> front_back.weight('aaa')
	2.0
	>>> front_back.weight('idyto')
	0.0
	'''
	classes = nat_class_set.classes
	class1 = frozenset(class1)
	class2 = frozenset(class2)
	#states: '#': no part of pattern processed
	#		 '1': character from class1 processed; looking for character from class2
	start = '#'
	stops = {'#':-1, '1':-1}
	arcs = {}
	arcs['#'] = {}
	arcs['#'][class1] = ('1', -1)
	arcs['#']['_other'] = ('#', -1)
	arcs['1'] = {}
	intersect = frozenset(class1.intersection(class2))
	arcs['1'][intersect] = ('1', 0)
	cover_labels = get_covering_labels(class1.difference(intersect), classes)
	for label in cover_labels: #arcs for characters in class1 and not class2
		arcs['1'][label] = ('1', -1)
	cover_labels = get_covering_labels(class2.difference(intersect), classes)
	for label in cover_labels: #arcs for characters in class2 and not class1
		arcs['1'][label] = ('#', 0)
	labels = get_covering_labels(tier.difference(class1).difference(class2), classes)
	for label in labels: #arcs
		arcs['1'][label] = ('#', -1)
	arcs['1']['_other'] = ('1', -1) #arc for all letters not on this tier
	return NaturalClassWFSA( nat_class_set, start, stops, arcs, [param] )

def class_trigram( nat_class_set, class1, class2, class3, tier = None, param = 1.0 ):
	'''
	>>> consonant = frozenset('ptbdszmn')
	>>> obst = frozenset('ptbdsz')
	>>> voice = frozenset('bdzmn')
	>>> vowel = frozenset('aeiou')
	>>> vless = frozenset('pts')
	>>> alphabet = set('ptbdszmnaeiou')
	>>> classes = NaturalClassSet(alphabet, set([consonant, obst, voice, vowel]))
	>>> obst_voice_C = class_trigram( classes, obst, voice, consonant )
	>>> obst_voice_C.weight('pbns')
	2.0
	>>> obst_voice_C.weight('tmnzabb')
	1.0
	>>> obst_voice_C.weight('tmsbtnm')
	3,0
	>>> obst_voice_vless = class_trigram( classes, obst, voice, vless, consonant )
	>>> obst_voice_vless.weight('pabos')
	1.0
	>>> obst_voice_vless.weight('pbasbazt')
	2.0
	'''
	classes = nat_class_set.classes
	class1 = frozenset(class1)
	class2 = frozenset(class2)
	class3 = frozenset(class3)
	class1and2 = class1.intersect(class2)
	class2and3 = class2.intersect(class3)
	class1and3 = class1.intersect(class3)
	#states: 	'#': no part of pattern processed
	#		 	'1': character from class1 processed; looking for character from class2
	#			'1&2': a seequence with one label from class1 and a second label in
	#					both class1 and class2 processed. 
	#			'2': a sequence with one label from class1 and one from class2 and not
	#					class1 has been processed.
	start = '#'
	stops = {'#':-1, '1':-1, '2':-1}
	arcs = {}
	arcs['#'] = {}
	arcs['#'][class1] = ('1', -1)
	arcs['#']['_other'] = ('#', -1)
	
	if len(class2and3) > 0 and len(class1and2) > 0:
		stops['1&2'] = -1
		arcs['1&2'] = {}
		arcs['1&2'][class2and3] = ('1&2', 0)
		labels = get_covering_labels( class1and3.difference(class2), classes )
		for label in labels:
			arcs['1

	arcs['1'] = {}

if __name__ == '__main__':
	import doctest
	doctest.testmod()
