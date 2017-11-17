import sys
sys.path.append("..")

from functools import reduce
from itertools import combinations,permutations
import ast

def get_class(polynomial, class_list, n=5):
	"""
	Gets the polynomial class this polynomial is equivalent to, or None if
	no valid class is found
	"""
	possible_permutations = list(permutations(range(n)))
	polynomial_class = [booleantools.apply_perm_to_polynomial(perm, polynomial) for perm in possible_permutations]
	for poly in class_list:
		for poss_class in polynomial_class:
			if booleantools.is_equal_polynomials(poly, poss_class):
				return poly

	return None
	
def get_poly_from_line(line):
	junk, arr = line.split(":")
	arr = arr.strip() # This string is a dirty one
	return ast.literal_eval(arr)

def get_res(poly):
	last_level = 0
	last_resilient = True
	while last_resilient:
		last_resilient = booleantools.is_k_resilient(booleantools.polynomial_to_rule_table(poly, 5), last_level + 1)
		if last_resilient:
			last_level += 1
	return last_level


def main():
	polynomial_classes = []
	cardinality = []
	with open("2res", "r") as file:
		for line in file:
			polynomial = get_poly_from_line(line)
			poly_class = get_class(polynomial, polynomial_classes)
			if poly_class is not None:
				index = polynomial_classes.index(poly_class)
				cardinality[index] += 1
			else:
				polynomial_classes.append(polynomial)
				cardinality.append(1)
			
	print("Analysis of results finished:")
	for polynomial_class, card in zip(polynomial_classes, cardinality):
		k = get_res(polynomial_class)
		print("Polynomial class", polynomial_class, "has", card, "members and is", k, "resilient")
	
	
	
	
if __name__ == "__main__":
	main()
