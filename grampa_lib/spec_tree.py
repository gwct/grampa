import sys
import os
import grampa_lib.reconcore as RC
import grampa_lib.recontree as RT

#############################################################################

def readSpecTree(globs, label_opt=False):
	if globs['st-input-type'] == "file":
		globs['st-str'] = open(globs['spec-tree-input'], "r").read().strip();
	# If the input string is a filename, read the file. Otherwise, just try it as a newick string.

	spec_tree = RT.remBranchLength(globs['st-str']);
	tips = spec_tree.replace("(","").replace(")","").replace(";","").split(",");
	# Remove the branch lengths from the tree string and get the tip labels.

	if any(tip.isdigit() for tip in tips):
		RC.errorOut("ST1", "Tip labels cannot be purely numbers. Please add another character.", globs);

	if globs['mul-input-flag'] and any(tips.count(tip) not in [1,2] for tip in tips):
		RC.errorOut("ST2", "You have entered a tree type (-t) of 'm', species in your tree should appear exactly once or twice.", globs);
	elif any(tips.count(tip) > 1 for tip in tips):
		RC.errorOut("ST3", "You have entered a tree type (-t) of 's' but there are labels in your tree that appear more than once!", globs);
	# Some error checking based on the tip labels in the tree.

	if globs['mul-input-flag']:
		hybrid_spec = list(set([tip for tip in tips if tips.count(tip) != 1]));
		for h in hybrid_spec:
			spec_tree = spec_tree.replace(h, h+"*", 1);
	# If the user entered a MUL-tree, some internal re-labeling must be done to those labels that appear twice.

	try:
		globs['st'], globs['parsed-st-str'] = RT.treeParse(spec_tree);
		# Parsing of the species tree.
	except:
		RC.errorOut("ST4", "Error reading species tree!");
	# Reading the species tree file.

	if label_opt:
		print();
		print("# The input species tree with internal nodes labeled:");
		print(globs['parsed-st-str'] + "\n");
		sys.exit(0);
	# The output if --labeltree is set.

	globs['tips'] = tips;

	return globs;

#############################################################################

def hInParse(globs):
	if globs['mul-input-flag']:
		mul_copy_clade = [n for n in sinfo if sinfo[n][2] == 'tip' and '*' in n];
		mul_hybrid_clade = [n.replace("*","") for n in mul_copy_clade];
		# Read the nodes with the added * to get the hybrid clades

		mul_hybrid_node, mul_hybrid_mono = RT.LCA(mul_hybrid_clade, sinfo);
		mul_copy_node, mul_copy_mono = RT.LCA(mul_copy_clade, sinfo);
		# Get the LCA for both clades

		if not mul_hybrid_mono or not mul_copy_mono:
			RC.errorOut("ST5", "All hybrid clades specified in your MUL-tree must be monophyletic! Hybrid clade identified as: " + ",".join(mul_copy_clade), globs);
		# Make sure the clades are monophyletic

		#hybrid_clades, hybrid_nodes, copy_clades, copy_nodes = [mul_hybrid_clade], [mul_hybrid_node], [mul_copy_clade], [mul_copy_node];
		globs['h1-clades'], globs['h1-nodes'], globs['h2-clades'], globs['h2-nodes'] = [mul_hybrid_clade], [mul_hybrid_node], [mul_copy_clade], [mul_copy_node];
		# Store the clades as variables
	# If the input tree is a MUL-tree, we have to determine what the hybrid clades and nodes are from the input tree.

	else:
		globs['h1-clades'], globs['h1-nodes'] = getHClades("h1", globs);
		globs['h2-clades'], globs['h2-nodes'] = getHClades("h2", globs);
	# If the input tree is singly-labeled, use the input info from -h1 and -h2 to get the hybrid clades and nodes.

	return globs;
# Parses the input h nodes.

#############################################################################

def getHClades(h_type, globs):
# This function takes a list of lists of -h1 or -h2 inputs and determines if they are clades or node labels. It then retrieves
# the complete lists of hybrid clades and nodes.
	h_list = globs[h_type + "-input"];

	if h_list:
		if " " in h_list:
			h_clades = h_list.split(" ");
			h_clades = list(map(set, [tmp_h.split(",") for tmp_h in h_clades]));
		else:
			h_clades = list(map(set, [h_list.split(",")]));
		# Split up the input info. If there is a space, multiple nodes/clades have been specified.


		missing_nodes = [];
		for hybrid_list in h_clades:
			for h in hybrid_list:
				if h.isdigit():
					h = "<" + h + ">";
				if h not in globs['st']:
					missing_nodes.append(h);

		if missing_nodes:
			RC.errorOut("ST6", "Not all -" + h_type + " species are present in your species tree: " + ", ".join(missing_nodes), globs);

		#if not all(h in sinfo for hybrid_list in h_clades for h in hybrid_list if not h.isdigit()):
		#	RC.errorOut("ST6", "Not all -" + h_type + " species are present in your species tree!");
		#if not all("<" + h + ">" in sinfo for hybrid_list in h_clades for h in hybrid_list if h.isdigit()):
		#	RC.errorOut("ST7", "Not all -" + h_type + " nodes are present in your species tree!");
		# Some error checking to make sure everything the user input is actually in the tree.

		h_nodes = [];
		for hybrid_clade in h_clades:
			hybrid_clade = list(hybrid_clade);
			if hybrid_clade[0].isdigit():
				h_node = "<" + hybrid_clade[0] + ">";
			# If the input was an internal node, add it to the node list here.
			else:
				h_node, h_mono = RT.LCA(hybrid_clade, globs['st']);
				if not h_mono:
					RC.errorOut("ST7", "All hybrid clades specified h1 and h2 must be monophyletic!", globs);
			# If the input was a clade, retrieve the ancestral node and check if it is monophyletic here.
			if h_node not in h_nodes:
				h_nodes.append(h_node);
			# Add the hybrid node to the nodes list.
	# If the user input anything as -h1 or -h2 this parses it.

	else:
		h_nodes = list(globs['st'].keys());
		h_clades = [RT.getClade(node, globs['st']) for node in h_nodes];
	# If the user did not specify -h1 or -h2, this adds all possible nodes to the list.

	return h_clades, h_nodes;

#############################################################################

