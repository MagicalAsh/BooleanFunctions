from enum import Enum as _Enum
import math as _math

class PrimeField:
    """
    Represents a finite field of prime order.
    """
    class _FieldElement:
        def __init__(self, value, field):
            self.value = value
            self.field = field

        def __add__(a, b):
            if a.field == b.field:
                return a.field.get((a.value+b.value) % a.field.order)

        def __mul__(a, b):
            if a.field == b.field:
                return a.field.get((a.value*b.value) % a.field.order)

        def __repr__(self):
            return "<%d>" % self.value
        
        def __str__(self):
            return "[%d]" % self.value

    def __init__(self, order):
        if not _isprime(order):
            raise ValueError("order is not prime")

        self.order = order
        self.elements = {ele : PrimeField._FieldElement(ele, self) for ele in range(0, order)}
        
    def get(self, element):
        if element > self.order or element < 0:
            raise ValueError("Element must be greater than zero and less than order")
        return self.elements[element]
    
    def value_of(self, element):
        if element not in self:
            raise ValueError("Element is not of this field")
        
        return element.value

    def __eq__(a, b):
        return True if a.order == b.order else False
    
    def __contains__(self, elem):
        if not isinstance(elem, PrimeField._FieldElement):
            return False

        if elem.field == self:
            return True

        return False

class GF4(_Enum):
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

def _isprime(n):
    if n < 2:
        return False
    elif n == 2 or n == 3:
        return True
    elif n % 2 == 0 or n % 3 == 0:
        return False
    
    for i in range(5, int(_math.sqrt(n))+1, 6):
        if n % i == 0 or n % (i+2) == 0:
            return False

    return True

