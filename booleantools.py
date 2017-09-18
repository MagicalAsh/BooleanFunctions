from functools import reduce
from itertools import combinations,permutations
from math import log

class BooleanFunction:
	def __init__(self, listform, n):
		self.n = n
		self.listform = listform
		self.update_rule_table()
	
	def hamming_weight(self):
		return sum(self.tableform)
		
	def hamming_distance(self, other):
		if hasattr(other, "__getitem__"): #If other is a list
			return [self.hamming_distance(f) for f in other]
		else: 
			u = self.tableform
			v = other.tableform
			s = sum([delta(u[k],v[k])for k in range(len(u))])
			return s
		
	def walsh_transform(self):
		f = self.tableform
		nbits = int(log(len(f),2))
		vecs = [dec_to_bin(x,nbits) for x in range(len(f))]
		def Sf(w):
		  return sum([(-1)**(f[x]^dot_product(dec_to_bin(x,nbits),w)) for x in range(0,len(f))])
		Sflist = [Sf(vec) for vec in vecs]
		return Sflist

	def walsh_spectrum(self):
		f = self.tableform
		walsh_transform_f = self.walsh_transform()
		spec = max([abs(v) for v in walsh_transform_f])
		return spec

	def is_balanced(self):
		f = self.tableform
		return sum(f) == len(f)/2

	def is_correlation_immune(self,k=1):
		f = self.tableform
		walsh_transform_f = self.walsh_transform()
		nbits = int(log(len(f),2))
		vectors_to_test = [bin_to_dec(vec) for vec in weight_k_vectors(k,nbits)]
		walsh_transform_at_weight_k = [walsh_transform_f[vec] for vec in vectors_to_test]
		return walsh_transform_at_weight_k == [0]*nbits

	def is_k_resilient(self,k=1):
		f = self.tableform
		return self.is_balanced() and self.is_correlation_immune(k)
	
	#nonlinearity using Walsh transform
	def nonlinearity(self,n):
		f = self.tableform
		return 2**(n-1) - 0.5*self.walsh_spectrum()
	
	def evaluate_polynomial_at_point(self, x):
		f = self.listform
		value = 0
		for monomial in f:
			monomial_eval = product([x[i] for i in monomial])
			value += monomial_eval
		if [] in f:
			 value += 1
		value = value%2
		return value
		
	def apply_permutation(self, perm):
		def apply_perm_to_monomial(perm,monomial):
			out = [perm[i] for i in monomial]
			return out
	
		out = [apply_perm_to_monomial(perm, i) for i in self.listform]
		return BooleanFunction(out, self.n)
		
	def __eq__(self,poly2):
		return frozenset([frozenset(i) for i in self.listform]) == frozenset([frozenset(j) for j in poly2.listform])

	def __add__(self, other):
		newPoly = BooleanFunction(self.listform + other.listform)
		return newPoly
		
	def __mult__(self, other):
		if isinstance(other, int ): #I do this explicitly because it's kinda weird otherwise
			other = BooleanFunction(self.listform)
			for monomial in other.listform:
				monomial.append(other)
			other.n += 1
			other.update_rule_table()
		else:
			return None
			
	def __str__(self):
		return str(self.listform)
	
	def __repr__(self):
		return "BooleanFunction(%s)" % (str(self))
	
	def update_rule_table(self):
		rule_table_length = 2**self.n
		rule_table = [0]*rule_table_length
		for k in range(rule_table_length):
			point_to_evaluate = dec_to_bin(k, self.n)
			rule_table[k] = self.evaluate_polynomial_at_point(point_to_evaluate)
		self.tableform = rule_table
		
		
		
		
		
		
##A basic binary-to-decimal converter.
##Could obviously be optimized to reduce exponentiation.

def bin_to_dec(num):
	return sum([num[i]*2**i for i in range(len(num))])

#A basic decimal-to-binary converter.
#We need nbits in case padded 0's are needed at the front. 

def dec_to_bin(num,nbits):
	new_num = num
	bin = []
	for j in range(nbits):
		#create the appropriate power of 2 for the current step
		current_bin_mark = 2**(nbits-1-j)
		#check to see if you can subtract this power of 2; if so, then subtract it and append 1
		if (new_num >= current_bin_mark):
			bin.append(1)
			new_num = new_num - current_bin_mark
	  #if you can't subtract, append 0
		else:
			bin.append(0)
	return bin

#delta = lambda x,y: x==y # NOTE: Boolean values are actually a subclass of integers, so True*3 == 3
def delta(x,y):
	return x != y

def hausdorff_distance_point(a,B):
	"""
	Calculates the minimum distance between function a and the functions in the set B.
	"""
	return min([a.hamming_distance(b) for b in B])

def hausdorff_semidistance_set(A,B):
	return max([hausdorff_distance_point(a,B) for a in A])

def hausdorff_distance_sets(X,Y):
	"""
	Calculates the Hausdorff distance between two sets of boolean functions
	"""
	HD1 = hausdorff_semidistance_set(X,Y)
	HD2 = hausdorff_semidistance_set(Y,X)
	return max([HD1,HD2])
	 
def dot_product(u,v):
	"""
	Basic mod 2 dot product.
	"""
	s = sum(u[k]*v[k] for k in range(len(u)))
	return s%2

#Generate vectors with weight k
def weight_k_vectors(k,nbits):
	"""
	Generates all vectors with weight k.
	"""
	nums = range(nbits)
	vector_set_to_return = []
	k_combinations = [list(x) for x in combinations(nums,k)]
	for j in k_combinations:
		vec_to_add = [int(y in j) for y in range(nbits)]
		vector_set_to_return.append(vec_to_add)
	return vector_set_to_return
	
def product(x):
	return reduce((lambda y,z : y*z), x)
#Return the powerset of a list

def powerset(iterable):
	s = list(iterable)
	set_to_return = []
	for r in range(len(s)+1):
		sr_combinations = combinations(s,r)
		for item in sr_combinations:
			set_to_return.append(list(item))
	return set_to_return

#Return the powerset of a list, excluding the empty set
def nonemptypowerset(iterable):
	s = list(iterable)
	set_to_return = []
	for r in range(1,len(s)+1):
		sr_combinations = combinations(s,r)
		for item in sr_combinations:
			set_to_return.append(list(item))
	return set_to_return	
	
def duplicate_free_list_polynomials(list_of_polys):
	outlist = [] 
	for poly in list_of_polys:
		if True not in [poly == poly_in_out for poly_in_out in outlist]:
			outlist.append(poly)
	return outlist

	#todo: update me
def perms_orbit_polynomial(permset,polynomial):
	return duplicate_free_list_polynomials([polynomial.apply_permutation(i) for i in permset])

def perms_equiv_classes_polynomial_list(permset,polynomial_list):
	return [perms_orbit_polynomial(permset,polynomial) for polynomial in polynomial_list]
				
	#todo: fix to use clone
def siegenthaler_combination(f1,f2,new_var):
	f1_times_new_var = f1 * new_var
	f2_times_one = f2
	f2_times_new_var = f2 * new_var
	return f1_times_new_var + f2_times_one + f2_times_new_var

def generate_all_seigenthaler_combinations(func_list,new_var):
	all_siegenthaler_combinations = []
	for f1 in func_list:
		for f2 in func_list:
			f1f2siegenthalercombination = siegenthaler_combination(f1,f2,new_var)
			all_sigenthaler_combinations.append(f1f2siegenthalercombination)
	return all_seigenthaler_combinations
			
def test_function_in_equiv_classes(func,eq_classes):
	for eq_class in eq_classes:
		func_in_eq_class = False
		func_in_eq_class = [func == x for x in eq_class]
		return func_in_eq_class
		
def min_nonzero_dist(poly1, classA):
	dists = [poly1.hamming_distance(f) for f in classA]
	min_nonzero = float("inf")
	for dist in dists:
		if dist != 0 and dist < min_nonzero:
			min_nonzero = dist
			
	return dist


