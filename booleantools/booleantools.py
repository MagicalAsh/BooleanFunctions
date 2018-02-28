import copy
from enum import Enum
from functools import reduce
from itertools import combinations,permutations
from math import log
from copy import deepcopy


def Sym(n):
    """
    Creates a set containing all permutations in the symmetric group $S_n$.

    Returns:
        list: A set containing every permutation in $S_n$, in one-line notation.
    """
    return list(permutations([i for i in range(n)]))


class GF4(Enum):
    """
    Represents the Galois Field GF(4) with proper addition and multiplication.
    """
    ZERO = (0, 0)
    ONE = (0, 1)
    X = (1, 0)
    X_PLUS_ONE = (1, 1)
    def __add__(a, b):
        return GF4(((a.value[0] + b.value[0])%2, (a.value[1] + b.value[1])%2))

    def __mul__(a, b):
        if a == GF4.ZERO or b == GF4.ZERO:
            return GF4.ZERO
        else:
            z_mod_3 = {GF4.ONE: 0, GF4.X: 1, GF4.X_PLUS_ONE: 2}
            flipped_z3 = {0: GF4.ONE, 1: GF4.X, 2: GF4.X_PLUS_ONE}
            a_val = z_mod_3[a]
            b_val = z_mod_3[b]
            final = flipped_z3[(a_val + b_val) % 3]
            return final

class GF2(Enum):
    """
    Represents the Finite Field $F_2$ with proper addition and multiplication.
    """
    ZERO = 0
    ONE = 1
    def __add__(a, b):
        return GF2((a.value + b.value) % 2)

    def __mul__(a, b):
        return GF2(a.value * b.value)

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
        for pos, val in enumerate(args): #Simplification for inputs
            if val == 1:
                args[pos] = self.field.ONE
            elif val == 0:
                args[pos] = self.field.ZERO

        value = self.field.ZERO
        for monomial in self.listform:
            if monomial not in self.field:
                prod = self.field.ONE
                for var in monomial:
                    if var not in self.field:
                        prod *= args[var]
                    else:
                        prod *= var
                value += prod
            else:
                value += monomial
        return value
        
    
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
        return FieldFunction(a.listform + b.listform, max(a.n, b.n), a.field)

    def __mul__(a, b):
        if a.field != b.field:
            raise ValueError("Multiplicands are not from the same field.")

        out = []
        for monomial in a.listform:
            out += [monomial + monomial_b for monomial_b in b.listform]
        
        return FieldFunction(out, max(a.n, b.n), a.field)

    
    def tex_str(self, math_mode=False):
        """
        Creates a TeX String from this BooleanFunction.

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
        return "FieldFunction(%s, %s, %s)" % (str(self.listform), str(self.n), str(self.field))

    def __reduce(self):
        for i in range(len(self.listform)):
            self.listform[i] = [val for val in self.listform[i] if val != self.field.ONE]


class BooleanFunction(FieldFunction):
    """
    This class represents a boolean function ($\\mathbb{F}_2^n \\rightarrow \\mathbb{F}_2) and implements a large amount of useful functions.
    """
    def __init__(self, listform, n):
        """
        Creates a boolean function on n variables.
        
        Attributes:
            listform (list): A list of the monomials this polynomial contains. Ex. \[x_1 \oplus x_2x_3\] is [[0], [1, 2]].
            n (int): The number of variables, where n - 1 is the highest term in the list form.
        """
        super().__init__(listform, n, GF2)
        copyList = []

        #This is done for space efficiency
        for i in listform:
            if i not in copyList:
                copyList.append(i)
            else:
                copyList.remove(i)

        self.listform = copyList
        self.update_rule_table()
    
    def hamming_weight(self):
        """
            Returns the Hamming Weight of this function.
            
            Returns:
                int: The hamming weight of this function.
        """
        return sum(_GF2_to_ints(self.tableform))
        
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
            u = _GF2_to_ints(self.tableform)
            v = _GF2_to_ints(other.tableform)
            s = sum([_delta(u[k],v[k]) for k in range(len(u))])
            return s
        
    def walsh_transform(self):
        """
        Performs a Walsh transform on this function.

        Returns:
            list: A list containing the walsh transform of this function.
        """
        f = _GF2_to_ints(self.tableform)
        nbits = self.n
        vecs = [(_dec_to_bin(x,nbits), x) for x in range(len(f))]
        def Sf(w):
            return sum([(-1)**(f[x]^_dot_product(vec,w)) for vec,x in vecs])
        Sflist = [Sf(vec) for vec,x in vecs]
        return Sflist

    def __walsh_alt(self):
        f = _GF2_to_ints(self.tableform)
        vecs = [(_dec_to_bin(x, self.n), x) for x in range(len(f))]
        def walsh_t(w):
            return 2**self.n - 2*sum(f[i] ^ _dot_product(vec, w) for vec,i in vecs)
        return [walsh_t(vec) for vec,x in vecs]

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
        f = _GF2_to_ints(self.tableform)
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
            
        # See:
            #is_correlation_immune
            #is_balanced
        """
        return self.is_balanced() and self.is_correlation_immune(k=k)
    
    def is_affine(self):
        return True if self.nonlinearity() == 0 else False
    
    def get_class(self, perms=None):
        perms = perms if perms is not None else Sym(self.n)
        return perms_orbit_polynomial(perms, self)
     
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
    
    def tex_str(self, math_mode=False):
        """
        Creates a TeX String from this BooleanFunction.

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
    
    def apply_permutation(self, perm):
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
                    thisRow.append(self.func([a[0]] + j[:]))
                else:
                    thisRow.append(None)
            self.adjMatr.append(thisRow)

class CellularAutomata:
    """
    A class representing a basic cellular automata using a boolean function as a rule.
    """

    def __init__(self, func, initialState):
        self.func = func
        self.state = []
        for cell in initialState:
            if cell == 1:
                self.state.append(func.field.ONE)
            elif cell == 0:
                self.state.append(func.field.ZERO)
            else:
                self.state.append(cell)

        self.time = []

    def update(self):
        """
        Updates this Cellular Automata to the next state, and appends the old state to the time
        attribute.
        """
        self.time.append(copy.deepcopy(self.state))
        out = []
        for i in range(0, len(self.state)):
            start = i - (self.func.n//2)
            end = (i + (self.func.n//2)) % len(self.state) + 1
            out.append(self.func(*_circle_slice(self.state, start, end)))

        self.state = out
        
         
    def pretty_print(self, colorMap=None):
        """
        Prints the state over time in a format that is easily visible in a Bash command line.
        """
        colorMap = colorMap if colorMap is not None else {self.func.field.ZERO: " ", self.func.field.ONE: "\033[7m \033[0m"} 
        for row in self.time:
            for cell in row:
                print(colorMap[cell], sep="", end="")

            print()

    def get_column(self, col):
        """
        Gets a column over time of the CA.

        Args:
            col (int): An integer representing the column to get

        Returns:
            list: The column over time, with index 0 being the initial state.
        """
        return [row[col] for row in self.time]

class TargetTransitions(Enum):
        """
            An Enum representing which states under which a target CA is allowed to change state.
        """
        TARGET_TO_TARGET = 0
        TARGET_TO_NON = 1
        NON_TO_TARGET = 2
        NON_TO_NON = 3
        
        def standard_revs():
            """
            Returns the standard revision values, allowing moving from Target state to Target State,
            Nontarget state to Nontarget state, and Nontarget State to Target State.

            Returns:
                TargetTransitions: The above states.
            """
            return TargetTransitions.TARGET_TO_TARGET, TargetTransitions.NON_TO_TARGET, TargetTransitions.NON_TO_NON

class TargetedCellularAutomata(CellularAutomata):
    def __init__(self, func, initialState, targetState, *args):
        super().__init__(func, initialState)
        self.targetState = targetState
        self.reversionValues = frozenset(args)

    def update(self):
        super().update()
        oldState = self.time[-1]
        for i in range(len(oldState)):
            change = _get_change_state(self.targetState[i], oldState[i], self.state[i])
            if change not in self.reversionValues:
                self.state[i] = oldState[i]

def getX(n):
    """
    Gets a list of all possible x_i in order, from 0 to n-1.
    """
    return [BooleanFunction([[i]], n) for i in range(0, n)]

def _gen_atomic(n, pos):
        prod = BooleanFunction([[GF2.ONE]], n)
        for position, val in enumerate(_dec_to_bin(pos, n)):
            if val == 1:
                f = BooleanFunction([[position]], n)
                prod *= f
            else:
                f = BooleanFunction([[position], [GF2.ONE]], n)
                prod *= f
        if prod.tableform[pos] != GF2.ONE:
           raise BaseException("Bad things happened!")
        return prod

def _GF2_to_ints(lst):
    return [1 if x == GF2.ONE else 0 for x in lst]


def generate_function(rule_no, n): 
    endFunc = BooleanFunction([], n)
    binary_list = _dec_to_bin(rule_no, 2**n)
    for pos, val in enumerate(binary_list[::-1]):
        if val == 1:
            endFunc += _gen_atomic(n, pos)

    return endFunc

def _get_change_state(target, curr, nex):
    is_curr_tar = (target == curr)
    is_nex_tar = (target == nex)

    if is_curr_tar and not is_nex_tar:
        return TargetTransitions.TARGET_TO_NON
    elif is_nex_tar and not is_curr_tar:
        return TargetTransitions.NON_TO_TARGET
    elif is_nex_tar and is_curr_tar:
        return TargetTransitions.TARGET_TO_TARGET
    else:
        return TargetTransitions.NON_TO_NON


def _circle_slice(lst, start, end):
    """
    Slices a list in a loop. Treats the list as a ring and slices from beginning to end,
    looping back to the beginning of the list if necessary.
    """
    if 0 <= start < end < len(lst):
        return lst[start:end]
    elif start < 0:
        return lst[start+len(lst):] + lst[:end]
    elif end >= len(lst):
        return lst[start:] + lst[0:end-len(lst)]
    elif start > end:
        return lst[start:] + lst[:end]
    elif start == end:
        return lst[start:] + lst[:end]

    print("SLICE FAILURE: ", lst, "FROM:", start, "END:", end)
    return []

def _bin_to_dec(num):
    
    """
    Converts a binary vector to a decimal number.
    """
    return sum([num[i]*2**i for i in range(len(num))])

#A basic decimal-to-binary converter.
#We need nbits in case padded 0's are needed at the front. 

def _dec_to_bin(num,nbits):
    """
    Creates a binary vector of length nbits from a number.
    """
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
    k_combinations = [list(x) for x in combinations(nums,k)]
    for j in k_combinations:
        vec_to_add = [int(y in j) for y in range(nbits)]
        vector_set_to_return.append(vec_to_add)
    return vector_set_to_return
   
def weight_k_or_less_vectors(k, nbits):
    output = []
    for i in range(0, k+1):
        output += weight_k_vectors(i, nbits)

    return output
    
def _product(x):
    """
    Calculates the product of all elements in $x$.
    
    # Arguments
        x - A list to find the product of.
        
    # Returns
        The product of the list.
    """
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

def perms_orbit_polynomial(permset,polynomial):
    """
        Orbits a polynomial using the given permutation set.
        
        Args:
            permset: A set of permutations to apply to the function

        Returns:
            A list of the polynomials created by the given orbits.
    """
    return duplicate_free_list_polynomials([polynomial.apply_permutation(i) for i in permset])

def perms_equiv_classes_polynomial_list(permset,polynomial_list):
    """
    Orbits a list of polynomials using the given permutation set.
    
    Returns:
        A list of lists of the polynomials created by the given orbits.
    """
    return [perms_orbit_polynomial(permset,polynomial) for polynomial in polynomial_list]
                

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

def reduce_to_classes(f_list, permset):
    """
    Reduces a list of functions to a list of function classes given a permutation set.
    """
    basic_polys = []
    flatten = lambda l: [item for sublist in l for item in sublist]
    for f in f_list:
        if f not in flatten([perms_orbit_polynomial(permset, basic) for basic in basic_polys]):
            basic_polys.append(f)
    return basic_polys
