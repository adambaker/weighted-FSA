from math import log

def immutable(self, *args):
	'''
	To make an object immutable, set
	self.__setattr__ = immutable
	'''
	raise TypeError('Cannot modify an immutable object.')

class Semiring(object):
	name = ''
		
	def __float__(self):
		return float(self._value)
		
	def __str__(self):
		return str(self._value)
	
	def __repr__(self):
		return repr(self._value) + ' in a '+self.name+'semiring.'
	
	def __cmp__(self, other):
		if self._value < other._value:
			return -1
		elif self._value > other._value:
			return 1
		else:
			return 0
	
	def __hash__(self):
		return self._value.__hash__()
	
class Probability(Semiring):
	name = 'probability '
	def __init__(self, value):
		if hasattr(value, '_value'):
			value = value._value
		Semiring.__setattr__(self, '_value', value)

	__setattr__ = immutable
	__delattr__ = immutable
		
	def to_tropical(self, base=2):
		return Tropical(-log(self._value, base))
	
	def to_probability(self, base=2):
		return self
		
	def __add__(self, other):
		return Probability(self._value + other._value)
	
	def __mul__(self, other):
		return Probability(self._value*other._value)
		
	def __div__(self, other):
		return Probability(self._value/other._value)
	
	@property
	def star(self):
		if self._value >= 0 and self._value < 1:
			return Probability(1/(1-self._value))
		else:
			return Probability(float('inf'))
	
Probability.zero = Probability(0.0)
Probability.one = Probability(1.0)

class Tropical(Semiring):
	name = 'tropical '
	def __init__(self, value):
		if hasattr(value, '_value'):
			value = value._value
		Semiring.__setattr__(self, '_value', value)
	
	__setattr__ = immutable
	__delattr__ = immutable
	
	def to_probability( self, base=2 ):
		return Probability( base**(-self._value) )
		
	def to_tropical( self, base=2 ):
		return self
	
	def __add__(self, other):
		return Tropical(min(self._value, other._value))
	
	def __mul__(self, other):
		return Tropical(self._value + other._value)
	
	def __div__(self, other):
		return Tropical(self._value - other._value)
		
	@property
	def star(self):
		return Tropical(0.0)

Tropical.zero = Tropical(float('inf'))
Tropical.one = Tropical(0.0)
