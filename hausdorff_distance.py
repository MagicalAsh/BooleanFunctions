import ast
import json
import booleantools as bt
from itertools import permutations

 
def get_poly_from_line(line):
	junk, arr = line.split(":")
	arr = arr.strip() # This string is a dirty one
	return ast.literal_eval(arr)

def main():
	# Import polynomials
	classes = []
	with open("2res_classes.txt", "r") as file:
		for line in file:
			classes.append(get_poly_from_line(line))
	
	# Hausdorff distance between classes
	perms = list(permutations(range(5)))
	output_dict = []
	
	for classA in classes:
		for classB in classes:
			if classA != classB:
				classA_set = bt.perms_orbit_polynomial(perms, classA)
				classA_table = [bt.polynomial_to_rule_table(a, 5) for a in classA_set]
				classB_set = bt.perms_orbit_polynomial(perms, classB)
				classB_table = [bt.polynomial_to_rule_table(b, 5) for b in classB_set]
				haus_dist = bt.hausdorff_distance_sets(classA_table, classB_table)
				output_dict.append({"first": classA, "second": classB, "distance": haus_dist})
	with open("hausdorff_distance.json", "w") as file:
		jsondata = json.dumps({"distances":output_dict}, indent=4)
		file.write(jsondata)
	
	
	
	# Diameter of classes 
	
	output_dict = []
	for classA in classes:
		max_dist = 0;
		classA_set = bt.perms_orbit_polynomial(perms, classA)
		classA_table = [bt.polynomial_to_rule_table(a, 5) for a in classA_set]
		diameter = max([bt.hausdorff_distance_point(poly1, classA_table) for poly1 in classA_table])
		diameter_data = {"class": classA, "diameter":diameter}
		output_dict.append(diameter_data)
		
	with open("diameter.json", "w") as file:
		jsondata = json.dumps({"diameters":output_dict}, indent=4)
		file.write(jsondata)
	
	
	
	
	
if __name__ == "__main__":
	main()