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
