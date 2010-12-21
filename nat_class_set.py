def features_to_classes( features ):
	'''Takes a dictionary from feature names to the set of phonemes with that feature.
	Returns a NaturalClassSet of natural classes induced from those features.'''
	alphabet = set([])
	for feature_phones in feaures.values():
		alphabet.update(feature_phones)
	return NaturalClassSet( alphabet, set(features.values()) )

class NaturalClassSet( object ):
	def __init__(self, alphabet, classes ):
		self.alphabet = set(alphabet)
		self.classes = self._fill_natural_classes(classes)
		for letter in self.alphabet:
			self.classes.add(frozenset([letter]))
		
		#potential arc labels (a natural class or '_other') are given 
		#probabilities in proportion to the number of phonemes in 
		#that label, with '_other' arcs treated as if it was the 
		#entire alphabet (thus it is the highest probability arc)
		#the weights are -log(pr(label))
		self.labels_len = reduce( lambda x,y: x+len(y), self.classes, 0.0 )
		self.labels_len += len(self.alphabet) #for the '_other' arc label
	
	def __iter__(self):
		return self.classes.__iter__()
	
	def __contains__(self, item):
		return item in self.classes
	
	def _fill_natural_classes(self, classes ):
		'''The intersection of any two natural classes should be a natural class.
		This function takes a set of sets (of phonemes) and intersects each set 
		with each other set, putting the result in a frozen set and adding it to the
		total set of natural classes. It does this until the set of classes stops
		growing, and returns the resulting set.
		
		This ensures that the set of natural classes is closed under set intersection.'''
		full_class_set = set([])
		start_size = len(classes)
		for c in classes:
			full_class_set.add(frozenset(c))
			for c2 in classes:
				full_class_set.add(frozenset(c.intersection(c2)))
		while start_size < len(full_class_set):
			start_size = len(full_class_set)
			classes = list(full_class_set)
			for c in classes:
				for c2 in classes:
					full_class_set.add(frozenset(c.intersection(c2)))
		full_class_set.remove(frozenset([]))
		return full_class_set
	
