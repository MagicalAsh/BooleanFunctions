from functools import reduce
from itertools import combinations,permutations
from math import log
from copy import deepcopy

class BooleanFunction:
    """
    This class represents a boolean function ($\mathbb{F}_2^n \\rightarrow \mathbb{F}_2$) and implements a large amount of useful functions.
    """
    def __init__(self, listform, n):
        """
        Creates a boolean function on $n$ variables.
        
        # Arguments
            listform -- A list of the monomials this polynomial contains. Ex. \[x_1 \oplus x_2x_3\] is [[0], [1, 2]].
            n -- The number of variables, where n - 1 is the highest term in the list form.
        """
        self.n = n
        self.listform = deepcopy(listform)
        self.__update_rule_table()
    
    def hamming_weight(self):
        """
            Returns the Hamming Weight of this function.
            
            # Returns
                The hamming weight of this function.
        """
        return sum(self.tableform)
        
    def hamming_distance(self, other):
        """
        Determines the hamming distance of a function or a list of functions.
        
        # Arguments
            other - function or list of functions to find distance to.
            
        # Returns
            A list of distances if #other is a list, or a float if #other is another function.
        """
        if not isinstance(self, BooleanFunction) and hasattr(other, "__getitem__"): #If other is a list
            return [self.hamming_distance(f) for f in other]
        else: 
            u = self.tableform
            v = other.tableform
            s = sum([_delta(u[k],v[k])for k in range(len(u))])
            return s
        
    def walsh_transform(self):
        """
        Performs a Walsh transform on this function.
        """
        f = self.tableform
        nbits = int(log(len(f),2))
        vecs = [_dec_to_bin(x,nbits) for x in range(len(f))]
        def Sf(w):
          return sum([(-1)**(f[x]^_dot_product(_dec_to_bin(x,nbits),w)) for x in range(0,len(f))])
        Sflist = [Sf(vec) for vec in vecs]
        return Sflist

    def walsh_spectrum(self):
        """
        Generates the Walsh spectrum for this function.
        """
        f = self.tableform
        walsh_transform_f = self.walsh_transform()
        spec = max([abs(v) for v in walsh_transform_f])
        return spec

    def is_balanced(self):
        """
        Determines whether this function is balanced or not.
        
        # Returns
            True if balanced, False otherwise.
        """
        f = self.tableform
        return sum(f) == len(f)/2

    def is_correlation_immune(self,k=1):
        """
        Determines if this function is $k$ correlation immune.
        
        # Arguments
            k - immunity level
        """
        f = self.tableform
        walsh_transform_f = self.walsh_transform()
        nbits = int(log(len(f),2))
        vectors_to_test = [_bin_to_dec(vec) for vec in weight_k_vectors(k,nbits)]
        walsh_transform_at_weight_k = [walsh_transform_f[vec] for vec in vectors_to_test]
        return walsh_transform_at_weight_k == [0]*len(vectors_to_test)

    def is_k_resilient(self,k=1):
        """
        Determines if this boolean function is $k$-resilient.
        
        # Arguments
            k - immunity level
            
        # See:
            #is_correlation_immune
            #is_balanced
        """
        f = self.tableform
        return self.is_balanced() and self.is_correlation_immune(k)
    
    def nonlinearity(self,n):
        """
        Gets the nonlinearity of this boolean function.
        
        # Returns
            Nonlinearity of this boolean function.
            
        """
        f = self.tableform
        return 2**(n-1) - 0.5*self.walsh_spectrum()
    
    def evaluate_polynomial_at_point(self, x):
        """
        Evaluates this function at $x$.
        
        # Arguments
            x - The value to evaluate at.
            
        # Returns
            An integer result to $f(x)$.
        """
        f = self.listform
        value = 0
        for monomial in f:
            monomial_eval = _product([x[i] for i in monomial])
            value += monomial_eval
        if [] in f:
             value += 1
        value = value%2
        return value
        
    def apply_permutation(self, perm):
        """
        Applies a permutation to this function.
        \\\[
            f^\sigma(x)
        \\\]
        
        where sigma is in one line notation.
        
        # Arguments
            perm - The permutation to apply, in one line notation. 
            
        # Returns
            A boolean function where the permutation was applied.
        """
        def apply_perm_to_monomial(perm,monomial):
            out = [perm[i] for i in monomial]
            return out
            
    
        out = [apply_perm_to_monomial(perm, i) for i in self.listform]
        return BooleanFunction(out, self.n)
        
    def __eq__(self,poly2):
        return self.tableform == poly2.tableform

    def __add__(self, other):
        newPoly = BooleanFunction(self.listform + other.listform, max(self.n, other.n))
        return newPoly
        
    def __mul__(self, other):
        if isinstance(other, int ): #I do this explicitly because it's kinda weird otherwise
            n = self.n
            if n == other:
                n += 1
            newf = BooleanFunction(self.listform, n)
            for monomial in newf.listform:
                monomial.append(other)
            newf.__update_rule_table()
            return newf
        else:
            return None
            
    def __str__(self):
        return str(self.listform)
    
    def __repr__(self):
        return "BooleanFunction(%s, %s)" % (str(self), str(self.n))
    
    def __update_rule_table(self):
        rule_table_length = 2**self.n
        rule_table = [0]*rule_table_length
        for k in range(rule_table_length):
            point_to_evaluate = _dec_to_bin(k, self.n)
            rule_table[k] = self.evaluate_polynomial_at_point(point_to_evaluate)
        self.tableform = rule_table
        
    def __hash__(self):
        return _bin_to_dec(self.tableform)
    
    def __call__(self, *args):
        return self.evaluate_polynomial_at_point(*args)
    
class FiniteStateMachine:
    """
    A class representing a finite state automata for a boolean function.
    """    

    def __init__(self, func):
        self.n = func.n
        self.func = func
        self.__generateGraph()

    def __generateGraph(self):
        self.adjMatr = []
        vertices = [_dec_to_bin(i, self.n - 1) for i in range(0, 2**(self.n - 1))]
        for i in range(0, (2**(self.n-1))):
            thisRow = []
            for j in vertices:
                a = _dec_to_bin(i, self.n - 1) 
                if a[1:] == j[:len(j) - 1]:
                    thisRow.append(self.func.evaluate_polynomial_at_point([a[0]] + j[:]))
                else:
                    thisRow.append(None)
            self.adjMatr.append(thisRow)

class CellularAutomata:
    
    def __init__(self, func, initialState):
        self.func = func
        self.state = initialState
        self.time = [initialState]

    def update(self):
        out = []
        for i in range(0, len(self.state)):
            start = i - (self.func.n//2)
            end = (i + (self.func.n//2)) % len(self.state) + 1
            out.append(self.func(_circle_slice(self.state, start, end)))

        self.state = out
        self.time.append(out)
    
    def get_column(self, col):
        return [row[col] for row in self.time]

##A basic binary-to-decimal converter.
##Could obviously be optimized to reduce exponentiation.

def _circle_slice(lst, start, end):
    if 0 <= start < end < len(lst):
        return lst[start:end]
    elif start < 0:
        return lst[start+len(lst):] + lst[:end]
    elif end >= len(lst):
        return lst[start:] + lst[0:end-len(lst)]
    elif start > end:
        return lst[start:] + lst[:end]
    return []

def _bin_to_dec(num):
    return sum([num[i]*2**i for i in range(len(num))])

#A basic decimal-to-binary converter.
#We need nbits in case padded 0's are needed at the front. 

def _dec_to_bin(num,nbits):
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

#__delta = lambda x,y: x==y # NOTE: Boolean values are actually a subclass of integers, so True*3 == 3
def _delta(x,y):

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
    Calculates the Hausdorff distance between two sets of boolean functions.
    """
    HD1 = hausdorff_semidistance_set(X,Y)
    HD2 = hausdorff_semidistance_set(Y,X)
    return max([HD1,HD2])
     
def _dot_product(u,v):
    """
    Basic mod 2 dot product.
    """
    s = sum(u[k]*v[k] for k in range(len(u)))
    return s%2

def weight_k_vectors(k,nbits):
    """
    Generates all vectors with hamming weight k.
    """
    nums = range(nbits)
    vector_set_to_return = []
    k_combinations = [list(x) for x in combinations(nums,k)]
    for j in k_combinations:
        vec_to_add = [int(y in j) for y in range(nbits)]
        vector_set_to_return.append(vec_to_add)
    return vector_set_to_return
    
def _product(x):
    """
    Calculates the __product of all elements in $x$.
    
    # Arguments
        x - A list to find the __product of.
        
    # Returns
        The __product of the list.
    """
    return reduce((lambda y,z : y*z), x)

    
def _powerset(iterable):
    """
    Generates the __powerset of an iterable.
    """
    s = list(iterable)
    set_to_return = []
    for r in range(len(s)+1):
        sr_combinations = combinations(s,r)
        for item in sr_combinations:
            set_to_return.append(list(item))
    return set_to_return

def _nonempty_powerset(iterable):
    """
    Generates the __powerset of an interable, excluding the empty set.
    """
    s = list(iterable)
    set_to_return = []
    for r in range(1,len(s)+1):
        sr_combinations = combinations(s,r)
        for item in sr_combinations:
            set_to_return.append(list(item))
    return set_to_return    
    
def duplicate_free_list_polynomials(list_of_polys):
    """
    Takes a list of boolean functions and generates a duplicate free list of polynomials.
    
    # Arguments
        list_of_polys - A list of polynomials.
        
    # Returns
        A duplicate free list of functions
    """
    outlist = [] 
    for poly in list_of_polys:
        if True not in [poly == poly_in_out for poly_in_out in outlist]:
            outlist.append(poly)
    return outlist

def perms_orbit_polynomial(permset,polynomial):
    """
        Orbits a polynomial using the given permutation set.
        
        # Returns
            A list of the polynomials created by the given orbits.
    """
    return duplicate_free_list_polynomials([polynomial.apply_permutation(i) for i in permset])

def perms_equiv_classes_polynomial_list(permset,polynomial_list):
    """
    Orbits a list of polynomials using the given permutation set.
    
    # Returns
        A list of lists of the polynomials created by the given orbits.
    """
    return [perms_orbit_polynomial(permset,polynomial) for polynomial in polynomial_list]
                

def siegenthaler_combination(f1,f2,new_var):
    """
    Generates a Siegenthaler Combination of two functions.
    
    # Arguments
        f1 - The first function
        f2 - The second function
        new_var - New variable for the combined function.
    
    # Returns
        The Siegenthaler combination of $f_1$ and $f_2$
        
    """
    f1_times_new_var = f1 * new_var
    f2_times_one = f2
    f2_times_new_var = f2 * new_var
    return f1_times_new_var + f2_times_one + f2_times_new_var

def generate_all_seigenthaler_combinations(func_list,new_var):
    """
    Generates all of the possible Siegenthaler combinations
    of the given functions.
    
    # Arguments
        func_list - A list of functions to perform the Siegenthaler combination function on.
        
    # Returns
        A list of all possible Siegenthaler combinations for the given functions.
    """
    all_siegenthaler_combinations = []
    for f1 in func_list:
        for f2 in func_list:
            f1f2siegenthalercombination = siegenthaler_combination(f1,f2,new_var)
            all_sigenthaler_combinations.append(f1f2siegenthalercombination)
    return all_seigenthaler_combinations
        
def min_nonzero_dist(poly1, classA):
    """
        Determines the minimum nonzero distance between a polynomial and its nearest neighbor.
        
        # Arguments
            poly1 - A boolean function
            classA - A class of boolean functions.
        
        # Returns
            The minimum nonzero distance between poly1 and every element of classA.
    """
    dists = [poly1.hamming_distance(f) for f in classA]
    min_nonzero = float("inf")
    for dist in dists:
        if dist != 0 and dist < min_nonzero:
            min_nonzero = dist
            
    return dist

def reduce_to_classes(f_list, permset):
    basic_polys = []
    flatten = lambda l: [item for sublist in l for item in sublist]
    for f in f_list:
        if f not in flatten([perms_orbit_polynomial(permset, basic) for basic in basic_polys]):
            basic_polys.append(f)
    return basic_polys
