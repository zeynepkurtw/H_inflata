from ete3 import Tree

# Define the Newick string
newick_string = "((kbiala:0.689519,carpe:0.841007)0.88208:0.129443,((wb:0.390814,muris:0.378621)0.98143:0.309776,((trepo:0.346537,HIN:0.31075)0.877437:0.190443,spiro:0.540709)0.821727:0.18251)0.88208:0.129443);"

# Create a Tree object from the Newick string
species_tree = Tree(newick_string)

# Print the species tree
print(species_tree)
