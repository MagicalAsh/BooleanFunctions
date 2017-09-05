import booleantools
import json
from itertools import permutations
from booleantools import BooleanFunction


def getNearestNeighbor(f, B):
	#There are multiple nearest neighbors
	"""
	Find the nearest neighbors in polynomial class B.
	"""
	nearest = []
	distance_to_nearest = float("inf")
	no_nearest = 0
	for g in B:
		distance_to_g = f.hamming_distance(g)
		if distance_to_g < distance_to_nearest:
			no_nearest = 1
			nearest = [g]
			distance_to_nearest = distance_to_g
		elif distance_to_g == distance_to_nearest :
			no_nearest += 1
			nearest.append(g)
			
	print("%d nearest neighbors, class B has %d members" % (no_nearest, len(B)))
	return nearest 

def main():
	perms = list(permutations(range(5)))
	
	polys = []
	with open("in/2res_classes.json", "r") as file:
		polys = [BooleanFunction(list_func, 5) for list_func in json.loads(file.read())["classes"]]
		
	print("This should have one element " + str(booleantools.duplicate_free_list_polynomials([polys[0],polys[0]])))
		
	# is rotation a distance preserving function?
	classes = booleantools.perms_equiv_classes_polynomial_list(perms, polys)
	
	#Uncomment this if you want to sit for a while
	
	for classA in classes:
		for classB in classes:						
			if classA != classB:
				for a in classA:
					nearest = getNearestNeighbor(a, classB)
					distance_before = a.hamming_distance(nearest)
					orbit_distances = []
					for perm in perms: #oh sweet jesus this is the worst algorithm I have ever written
						orbited_other = [f.apply_permutation(perm) for f in nearest]
						a_orbited = a.apply_permutation(perm)
						orbit_distances.extend(a_orbited.hamming_distance(orbited_other))
					if len(set(orbit_distances)) > 1:
						print(str(orbit_distances))
						print("FALSE")
	
	print("done")
	
	# Result: Nearest neighbors orbit together! Yay!
	
if __name__ == "__main__":
	main()