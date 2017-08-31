import ast
import json
import booleantools
from booleantools import BooleanFunction
from itertools import permutations

 
def get_poly_from_line(line):
	junk, arr = line.split(":")
	arr = arr.strip() # This string is a dirty one
	return ast.literal_eval(arr)

def main():
	# Import polynomials
	classes = []
	with open("in/2res_classes.json", "r") as file:
		classes = [BooleanFunction(list_func, 5) for list_func in json.loads(file.read())["classes"]]
	
	# Hausdorff distance between classes
	perms = list(permutations(range(5)))
	output_dict = []
	
	for classA in classes:
		for classB in classes:
			if classA != classB:
				classA_set = booleantools.perms_orbit_polynomial(perms, classA)
				classB_set = booleantools.perms_orbit_polynomial(perms, classB)
				haus_dist = booleantools.hausdorff_distance_sets(classA_set, classB_set)
				output_dict.append({"first": classA.listform, "second": classB.listform, "distance": haus_dist})
	with open("out/hausdorff_distance.json", "w") as file:
		jsondata = json.dumps({"distances":output_dict}, indent=4)
		file.write(jsondata)
		
	
	for classA in classes:
		for classB in classes:
			classA_set = booleantools.perms_orbit_polynomial(perms, classA)
			classB_set = booleantools.perms_orbit_polynomial(perms, classB)
			
			print("Classes %s and %s: " % (str(classA), str(classB)))
			
			for a in classA_set:
				# find a's nearest neighbor in b
				dist = booleantools.hausdorff_distance_point(a, classB_set)
				print(dist)
				
			
	
	
	# Diameter of classes 
	
	output_dict = []
	for classA in classes:
		max_dist = 0;
		classA_set = booleantools.perms_orbit_polynomial(perms, classA)
		diameter = max([booleantools.hausdorff_distance_point(poly1, classA_set) for poly1 in classA_set])
		diameter_data = {"class": classA.listform, "diameter":diameter}
		output_dict.append(diameter_data)
		
	with open("out/diameter.json", "w") as file:
		jsondata = json.dumps({"diameters":output_dict}, indent=4)
		file.write(jsondata)
	
	
if __name__ == "__main__":
	main()