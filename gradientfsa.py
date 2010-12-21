import wfsa as fsa
from optimize.gradient import Minimizer
import random
from math import log

class BoltzmannMinimizer( object ):
	'''BotlsmannMinimizer defines a params attribute and an objective method so
	that it may be used with optimize.gradient.Minimizer. The objective function
	is the cost in bits of encoding a random subset (chosen at initialization) of 
	the supplied corpus.'''
	
	def __init__(self, cond_bigrams, vowel_mi, corpus_list,
			min_cond = 0.01, min_char_encode = 0.001, num_words=5000):
		self.bigram_size = len(cond_bigrams.keys())
		self.bigram_keys = []
		self.vowel_keys = []
		self.params = []
		
		for bigram in cond_bigrams.keys():
			self.bigram_keys.append(bigram)
			self.params.append(cond_bigrams[bigram])
		
		for v_bigram in vowel_mi.keys():
			self.vowel_keys.append(v_bigram)
			self.params.append(vowel_mi[v_bigram])
		
		self.min_cond = min_cond
		self.min_char_encode = min_char_encode
		self.corpus = random.sample(corpus_list, num_words)
	
	def validate(self):
		'''Randomly searching the parameter space may result in an invalid parameter
		list in self.params. Validate checks each parameter value in self.params and
		changes each invalid value to the nearest valid parameter value.'''
		#ensure the conditional bigrams meet the minimum
		for i in range(self.bigram_size):
			if self.params[i] < self.min_cond:
				self.params[i] = self.min_cond
		
		
		for v_index, v_bigram in enumerate(self.vowel_keys):
			
			#constraint on the maximum value of this parameter
 			param_max = self.params[self.bigram_size+v_index] 
			
			for i, bigram in enumerate(self.bigram_keys):
				if v_bigram[1] == bigram[1]:
					if( self.params[i] - self.params[self.bigram_size+v_index]
							< self.min_char_encode ):
						param_max = min(param_max, self.params[i] - self.min_char_encode)
			
			self.params[self.bigram_size+v_index] = param_max
	
	def objective(self, params):
		self.validate()
		
		conditionals = {}
		for i in range(len(self.bigram_keys)):
			conditionals[self.bigram_keys[i]] = params[i]
			
		v_mi = {}
		for i in range(len(self.vowel_keys)):
			v_mi[self.vowel_keys[i]] = params[self.bigram_size+i]
		
		bigram_model = fsa.bigram_plog_wfsa( conditionals )
		
		alphabet = set([x[0] for x in self.bigram_keys])
		v_mi_model = fsa.v_mi_wfsa( alphabet, v_mi )
		v_mi_model = v_mi_model.intersect( bigram_model )
		v_mi_model = v_mi_model.mult_wfsa()
		v_mi_norm, converged = v_mi_model.norm_constant('#%#', 1e-6, 200)
		
		if not converged:
			return float('inf') #returning infinity ensures this solution isn't used.
			
		v_mi_len = 0
		for word in self.corpus:
			v_mi_len += -log(v_mi_model.weight(word)/v_mi_norm, 2)
		
		return v_mi_len
	
