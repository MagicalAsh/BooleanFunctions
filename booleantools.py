from functools import reduce
from itertools import combinations,permutations
from math import log
import pickle 


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

delta = lambda x,y: x==y # NOTE: Boolean values are actually a subclass of integers, so True*3 == 3

def hamming_distance(u,v):
    s = sum([delta(u[k],v[k])for k in range(len(u))])
    return s

def hamming_weight(u):
    return sum(u)

def hausdorff_distance_point(a,B):
    return min([hamming_distance(a,b) for b in B])

def hausdorff_semidistance_set(A,B):
    return max([hausdorff_distance_point(a,B) for a in A])

def hausdorff_distance_sets(X,Y):
    HD1 = hausdorff_semidistance_set(X,Y)
    HD2 = hausdorff_semidistance_set(Y,X)
    return max([HD1,HD2])
     
#A basic mod-two dot product
def dot_product(u,v):
    s = sum(u[k]*v[k] for k in range(len(u)))
    return s%2

#Generate vectors with weight k
def weight_k_vectors(k,nbits):
    nums = range(nbits)
    vector_set_to_return = []
    k_combinations = [list(x) for x in combinations(nums,k)]
    for j in k_combinations:
        vec_to_add = [int(y in j) for y in range(nbits)]
        vector_set_to_return.append(vec_to_add)
    return vector_set_to_return
    
#The guts of the Walsh transform
def walsh_transform(f):
     nbits = int(log(len(f),2))
     vecs = [dec_to_bin(x,nbits) for x in range(len(f))]
     def Sf(w):
          return sum([(-1)**(f[x]^dot_product(dec_to_bin(x,nbits),w)) for x in range(0,len(f))])
     Sflist = [Sf(vec) for vec in vecs]
     return Sflist

def walsh_spectrum(f):
    walsh_transform_f = walsh_transform(f)
    spec = max([abs(v) for v in walsh_transform_f])
    return spec

#Tests using Walsh transform
def is_balanced(f):
    return sum(f) == len(f)/2

def is_correlation_immune(f,k=1):
    walsh_transform_f = walsh_transform(f)
    nbits = int(log(len(f),2))
    vectors_to_test = [bin_to_dec(vec) for vec in weight_k_vectors(k,nbits)]
    walsh_transform_at_weight_k = [walsh_transform_f[vec] for vec in vectors_to_test]
    if walsh_transform_at_weight_k == [0]*nbits:
        return True
    else:
        return False

def is_k_resilient(f,k=1):
    if is_balanced(f):
        if is_correlation_immune(f,k):
            return True
        else:
            return False
    else:
        return False

#nonlinearity using Walsh transform
def nonlinearity(f,n):
    return 2**(n-1) - 0.5*walsh_spectrum(f)
    

#Small lambda function to multiply everything in a list together
product = lambda z: reduce((lambda x,y : x*y), z)

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
    
#Code to evaluate a function f at a point
def evaluate_polynomial_at_point(f,x):
    value = 0
    for monomial in f:
        monomial_eval = product([x[i] for i in monomial])
        value += monomial_eval
    if [] in f:
         value += 1
    value = value%2
    return value

#Convert from algebraic normal (polynomial) form to "rule table" form
def polynomial_to_rule_table(f, n = 4):
    rule_table_length = 2**n
    rule_table = [0]*rule_table_length
    for k in range(rule_table_length):
        point_to_evaluate = dec_to_bin(k,n)
        rule_table[k] = evaluate_polynomial_at_point(f,point_to_evaluate)
    return rule_table

def apply_perm_to_monomial(perm,monomial):
    out = [perm[i] for i in monomial]
    return out

def apply_perm_to_polynomial(perm,polynomial):
    out = [apply_perm_to_monomial(perm, i) for i in polynomial]
    return out
    
def is_equal_polynomials(poly1,poly2):
    return frozenset([frozenset(i) for i in poly1]) == frozenset([frozenset(j) for j in poly2])
    
def duplicate_free_list_polynomials(list_of_polys):
    outlist = [] 
    for poly in list_of_polys:
        if not(True in [is_equal_polynomials(poly,test) for test in outlist]):
            outlist.append(poly)
    return outlist

def perms_orbit_polynomial(permset,polynomial):
    return duplicate_free_list_polynomials([apply_perm_to_polynomial(i, polynomial) for i in permset])

def perms_equiv_classes_polynomial_list(permset,polynomial_list):
    return [perms_orbit_polynomial(permset,polynomial) for polynomial in polynomial_list]
                
def siegenthaler_combination(f1,f2,new_var):
    f1_times_new_var = [monomial.append(new_var) for monomial in f1]
    f2_times_one = list(f2)
    f2_times_new_var = [monomial.append(new_var) for monomial in f2]
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
        func_in_eq_class = [is_equal_polynomial(func,x) for x in eq_class]
        if func_in_eq_class:
            return True
        return func_in_eq_class


