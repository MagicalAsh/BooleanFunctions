import copy as _copy
import booleantools.fields as _fields
from itertools import combinations as _combs
from itertools import permutations as _perms

def Sym(n):
    """
    Creates a set containing all permutations in the symmetric group $S_n$.

    Returns:
        list: A set containing every permutation in $S_n$, in one-line notation.
    """
    return list(_perms([i for i in range(n)]))

GF2 = _fields.PrimeField(2)

class FieldFunction:
    """
    Represents a function over a Finite Field.
    """
    
    def __init__(self, listform, n, field):
        self.listform = listform
        self.n = n
        self.field = field
        self.__reduce()

    def __call__(self, *args):
        args = list(args)
        if len(args) == 1 and hasattr(args[0], '__getitem__'):
            args = args[0]
        elif len(args) == 1 and isinstance(args[0], int):
            args = _dec_to_base(args[0], self.n, self.field.order)

        for pos, val in enumerate(args): #Simplification for inputs
            args[pos] = self.field.get(val)
        
        value = self.field.get(0)
        for monomial in self.listform:
            if monomial not in self.field:
                prod = self.field.get(1)
                for var in monomial:
                    if var not in self.field:
                        prod *= args[var]
                    else:
                        prod *= var
                value += prod
            else:
                value += monomial
        
        return self.field.value_of(value)
        
    
    def apply_permutation(self, perm):
        """
        Applies a permutation to this function.
        \\\[
            f^\sigma(x)
        \\\]
        
        where sigma is in one line notation.
        
        Args:
            perm (list): The permutation to apply, in one line notation. 
            
        Returns:
            FieldFunction: A function where the permutation was applied.
        """
        def apply_perm_to_monomial(perm,monomial):
            out = []
            for var in monomial:
                if var in self.field:
                    out.append(var)
                else:
                    out.append(perm[var])
            return out
            
    
        out = [apply_perm_to_monomial(perm, i) for i in self.listform]
        return FieldFunction(out, self.n, self.field)

    def __add__(a, b):
        if a.field != b.field:
            raise ValueError("Summands from different fields.")
        if isinstance(b, _fields.PrimeField._FieldElement):
            return FieldFunction(a.listform + [[b]], a.n, a.field)
        else:
            return FieldFunction(a.listform + b.listform, max(a.n, b.n), a.field)

    def __mul__(a, b):
        if a.field != b.field:
            raise ValueError("Multiplicands are not from the same field.")
        
        if isinstance(b, _fields.PrimeField._FieldElement):
            out = [monomial + [b] for monomial in a.listform]
            return FieldFunction(out, a.n, a.field)
        else:
            out = []
            for monomial in a.listform:
                out += [monomial + monomial_b for monomial_b in b.listform]
        
            return FieldFunction(out, max(a.n, b.n), a.field)

    
    def tex_str(self, math_mode=False):
        """
        Creates a TeX String from this FieldFunction.

        Args:
            math_mode (bool, optional): Whether to return with surrounding '$'.

        Returns:
            str: A proper TeX String representing this function.
        """
        out = "" if not math_mode else "$"
        flag = False
        for monomial in self.listform:
            out += " \\oplus " if flag else ""
            for term in monomial:
                out += "x_{" + str(term) + "}"

            flag = True

        return out if not math_mode else out + "$"

    def __str__(self):
        return self.tex_str()

    def __repr__(self):
        return "FieldFunction(%s, %s, %s)" % (str(self.listform), str(self.n), \
                str(self.field))

    def __reduce(self):
        for i in range(len(self.listform)):
            self.listform[i] = [val for val in self.listform[i] if val != self.field.get(1)]


class BooleanFunction(FieldFunction):
    """
    This class represents a boolean function ($\\mathbb{F}_2^n \\rightarrow 
    \\mathbb{F}_2) and implements a large amount of useful functions.
    """
    def __init__(self, listform, n):
        """
        Creates a boolean function on n variables.
        
        Attributes:
            listform (list): A list of the monomials this polynomial contains. 
                             Ex. \[x_1 \oplus x_2x_3\] is [[0], [1, 2]].
            n (int): The number of variables, where n - 1 is the highest 
                     term in the list form.
        """
        super().__init__(listform, n, GF2)
        _copyList = []

        #This is done for space efficiency. Basically reduces coefficient mod 2
        for i in listform:
            if i not in _copyList:
                _copyList.append(i)
            else:
                _copyList.remove(i)

        self.listform = _copyList
        self.update_rule_table()
    
    def hamming_weight(self):
        """
            Returns the Hamming Weight of this function.
            
            Returns:
                int: The hamming weight of this function.
        """
        return sum(self.tableform)
        
    def hamming_distance(self, other):
        """
        Determines the hamming distance of a function or a list of functions.
        
        Args:
            other (BooleanFunction): function or list of functions to find distance to.
            
        Returns:
            int: A list of distances if #other is a list, or a float if #other is another function.
        """
        if hasattr(other, "__getitem__"): #If other is a list
            return [self.hamming_distance(f) for f in other]
        else: 
            u = self.tableform
            v = other.tableform
            s = sum([_delta(u[k],v[k]) for k in range(len(u))])
            return s
        
    def walsh_transform(self):
        """
        Performs a Walsh transform on this function.

        Returns:
            list: A list containing the walsh transform of this function.
        """
        f = self.tableform
        nbits = self.n
        vecs = [(_dec_to_bin(x,nbits), x) for x in range(len(f))]
        def Sf(w):
            return sum([(-1)**(f[x]^_dot_product(vec,w)) for vec,x in vecs])
        Sflist = [Sf(vec) for vec,x in vecs]
        return Sflist

    def walsh_spectrum(self):
        """
        Generates the Walsh spectrum for this function.

        Returns:
            float: The Walsh spectrum of this function.
        """
        f = self.tableform
        walsh_transform_f = self.walsh_transform()
        spec = max([abs(v) for v in walsh_transform_f])
        return spec

    def is_balanced(self):
        """
        Determines whether this function is balanced or not.
        
        # Returns
            bool: True if balanced, False otherwise.
        """
        f = self.tableform
        return sum(f) == len(f)/2

    def is_correlation_immune(self,k=1):
        """
        Determines if this function is k correlation immune.
        
        Args:
            k (int): immunity level
        """
        if k > self.n:
            raise BaseException("Correlation immunity level cannot be higher than the number of variables.")
        f = self.tableform
        walsh_transform = self.walsh_transform()
        nbits = self.n
        vectors_to_test = [_bin_to_dec(vec) for vec in weight_k_or_less_vectors(k,nbits)]
        walsh_transform_at_weight_k = [walsh_transform[vec] for vec in vectors_to_test]
        return walsh_transform_at_weight_k == [0]*len(walsh_transform_at_weight_k)

    def is_k_resilient(self,k=1):
        """
        Determines if this boolean function is k-resilient.
        
        Args:
            k (int): immunity level
        """
        return self.is_balanced() and self.is_correlation_immune(k=k)
    
    def is_affine(self):
        """
        Determines if this function is affine.

        Returns:
            True if this function is affine, false otherwise.
        """
        return True if self.nonlinearity() == 0 else False
    
    def get_orbit(self, perms=None):
        """
        Gets the orbit of this function under action of the symmetric group.

        Args:
            perms - default None. Uses this as a permutation set, otherwise the
                    full symmetric group on n symbols.

        Returns:
            A list containing all functions in the orbit of this function.
        """
        return orbit_polynomial(self, perms)
     
    def nonlinearity(self):
        """
        Gets the nonlinearity of this boolean function.
        
        Returns:
            int: Nonlinearity of this boolean function.
            
        """
        return int(2**(self.n-1) - 0.5*self.walsh_spectrum())
    
    def linear_structures(self):
        """
            Creates a set of values that exist as linear structures of this polynomial.

            Returns:
                set: Set of linear structures.
        """
        flatten = lambda l: [item for sublist in l for item in sublist]
        linear_structs = set(flatten(self.listform))
        for monomial in self.listform:
            if len(monomial) > 1:
                linear_structs -= set(monomial)

        return linear_structs
    
    def apply_permutation(self, perm):
        """
        Applies a permutation to the ordering of the variables to this function.

        Args:
            perm - The permutation to apply.

        Returns:
            The newly permuted function.
        """
        f = super().apply_permutation(perm)
        return BooleanFunction(f.listform, f.n)

    def __str__(self):
        return self.tex_str()
    
    def __add__(a, b):
        sum_f = FieldFunction.__add__(a, b)
        return BooleanFunction(sum_f.listform, sum_f.n)
    
    def __mul__(a, b):
        prod_f = FieldFunction.__mul__(a, b)
        return BooleanFunction(prod_f.listform, prod_f.n)

    def __eq__(self,poly2):
        return self.tableform == poly2.tableform

    def __repr__(self):
        return "BooleanFunction(%s, %s)" % (str(self.listform), str(self.n))
    
    def update_rule_table(self):
        rule_table_length = 2**self.n
        rule_table = [0]*rule_table_length
        for k in range(rule_table_length):
            point_to_evaluate = _dec_to_bin(k, self.n)
            rule_table[k] = self(*point_to_evaluate)
        self.tableform = rule_table
        
    def __hash__(self):
        return _bin_to_dec(self.tableform)
    

def getX(n, field=GF2):
    """
    Gets a list of all possible x_i in order, from 0 to n-1.
    """
    if field == GF2:
        return [BooleanFunction([[i]], i+1) for i in range(0, n)]
    else:
        return [FieldFunction([[i]], i+1, field) for i in range(0, n)]

def _gen_atomic(n, pos):
        prod = BooleanFunction([[GF2.get(1)]], n)
        for position, val in enumerate(_dec_to_bin(pos, n)):
            if val == 1:
                f = BooleanFunction([[position]], n)
                prod *= f
            else:
                f = BooleanFunction([[position], [GF2.get(1)]], n)
                prod *= f
        if prod.tableform[pos] != 1:
           raise BaseException("_gen_atomic failed! Please report on Github!")
        return prod

def _GF2_to_ints(lst):
    return [1 if x == GF2.get(1) else 0 for x in lst]


def generate_function(rule_no, n): 
    endFunc = BooleanFunction([], n)
    binary_list = _dec_to_bin(rule_no, 2**n)
    for pos, val in enumerate(binary_list[::-1]):
        if val == 1:
            endFunc += _gen_atomic(n, pos)

    return endFunc

def _bin_to_dec(num):
    
    """
    Converts a binary vector to a decimal number.
    """
    return sum([num[i]*2**i for i in range(len(num))])

def _dec_to_bin(num,nbits):
    """
    Creates a binary vector of length nbits from a number.
    """
    new_num = num
    bin = []
    for j in range(nbits):
        current_bin_mark = 2**(nbits-1-j)
        if (new_num >= current_bin_mark):
            bin.append(1)
            new_num = new_num - current_bin_mark
        else:
            bin.append(0)
    return bin

def _dec_to_base(num, nbits, base):
    """
    Creates a binary vector of length nbits from a number.
    """
    new_num = num
    bin = []
    for j in range(nbits):
        current_bin_mark = base**(nbits-1-j)
        if (new_num >= current_bin_mark):
            bin.append(1)
            new_num = new_num - current_bin_mark
        else:
            bin.append(0)
    return bin

def _delta(x,y):
    """
    Returns 1 if x and y differ, 0 otherwise.
    """
    return x != y

def _hausdorff_distance_point(a,B):
    """
    Calculates the minimum distance between function a and the functions in the set B.
    """
    return min([a.hamming_distance(b) for b in B])

def _hausdorff_semidistance_set(A,B):
    return max([_hausdorff_distance_point(a,B) for a in A])

def hausdorff_distance(X,Y):
    """
    Calculates the Hausdorff distance between two sets of boolean functions.
    """
    HD1 = _hausdorff_semidistance_set(X,Y)
    HD2 = _hausdorff_semidistance_set(Y,X)
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
    k_combinations = [list(x) for x in _combs(nums,k)]
    for j in k_combinations:
        vec_to_add = [int(y in j) for y in range(nbits)]
        vector_set_to_return.append(vec_to_add)
    return vector_set_to_return
   
def weight_k_or_less_vectors(k, nbits):
    """
    Generates all vectors of weight k on nbits bits.

    Args:
        k - weight
        nbits - the number of bits

    Returns:
        All vectors of weight k on nbits bits.
    """
    output = []
    for i in range(0, k+1):
        output += weight_k_vectors(i, nbits)

    return output
    
def _product(x):
    return reduce((lambda y,z : y*z), x)

def duplicate_free_list_polynomials(list_of_polys):
    """
    Takes a list of boolean functions and generates a duplicate free list of polynomials.
    
    # Arguments
        list_of_polys (BooleanFunction): A list of polynomials.
        
    # Returns
        list: A duplicate free list of functions
    """
    outlist = [] 
    for poly in list_of_polys:
        if True not in [poly == poly_in_out for poly_in_out in outlist]:
            outlist.append(poly)
    return outlist

def orbit_polynomial(polynomial, permset=None):
    """
        Orbits a polynomial using the given permutation set.
        
        Args:
            permset: A set of permutations to apply to the function

        Returns:
            A list of the polynomials created by the given orbits.
    """
    if permset is None:
        permset = Sym(polynomial.n)
    return duplicate_free_list_polynomials([polynomial.apply_permutation(i) for i in permset])

def orbit_polynomial_list(polynomial_list, permset=None):
    """
    Orbits a list of polynomials using the given permutation set.
    
    Returns:
        A list of lists of the polynomials created by the given orbits.
    """
    return [orbit_polynomial(polynomial, permset) for polynomial in polynomial_list]
                

def siegenthaler_combination(f1,f2,new_var):
    """
    Generates a Siegenthaler Combination of two functions.
    
    Args:
        f1 (BooleanFunction): The first function
        f2 (BooleanFunction): The second function
        new_var (int): New variable for the combined function.
    
    Returns:
        The Siegenthaler combination of $f_1$ and $f_2$
        
    """
    f1_times_new_var = f1 * new_var
    f2_times_one = f2
    f2_times_new_var = f2 * new_var
    return f1_times_new_var + f2_times_one + f2_times_new_var

def generate_all_siegenthaler_combinations(func_list,new_var):
    """
    Generates all of the possible Siegenthaler combinations
    of the given functions.
    
    Args:
        func_list - A list of functions to perform the Siegenthaler combination function on.
        
    Returns:
        A list of all possible Siegenthaler combinations for the given functions.
    """
    all_siegenthaler_combinations = []
    for f1 in func_list:
        for f2 in func_list:
            f1f2siegenthalercombination = siegenthaler_combination(f1,f2,new_var)
            all_siegenthaler_combinations.append(f1f2siegenthalercombination)
    return all_siegenthaler_combinations
        
def min_nonzero_dist(poly1, classA):
    """
    Determines the minimum nonzero distance between a polynomial and its nearest neighbor.
        
    Args:
        poly1 - A boolean function
        classA - A class of boolean functions.
        
    Returns:
        The minimum nonzero distance between poly1 and every element of classA.
    """
    dists = [poly1.hamming_distance(f) for f in classA]
    min_nonzero = float("inf")
    for dist in dists:
        if dist != 0 and dist < min_nonzero:
            min_nonzero = dist
            
    return dist

def reduce_to_orbits(f_list, permset):
    """
    Reduces a list of functions to a list of function classes given a permutation set.
    """
    basic_polys = []
    flatten = lambda l: [item for sublist in l for item in sublist]
    for f in f_list:
        if f not in flatten([orbit_polynomial(permset, basic) for basic in basic_polys]):
            basic_polys.append(f)
    return basic_polys
