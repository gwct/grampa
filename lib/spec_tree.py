import sys, os, reconcore as RC, recontree as RT, global_vars as globs

#############################################################################

def readSpecTree(spec_tree_input, starttime):
	if os.path.isfile(spec_tree_input):
		spec_tree = open(spec_tree_input, "r").read().replace("\n", "").replace("\r","");
	else:
		spec_tree = spec_tree_input;
	# If the input string is a filename, read the file. Otherwise, just try it as a newick string.
	hybrid_spec = "";

	spec_tree = RT.remBranchLength(spec_tree);
	tips = spec_tree.replace("(","").replace(")","").replace(";","").split(",");
	# Remove the branch lengths from the tree string and get the tip labels.

	if any(tip.isdigit() for tip in tips):
		RC.errorOut(6, "Tip labels cannot be purely numbers. Please add another character.");
	if globs.spec_type == 's' and any(tips.count(tip) > 1 for tip in tips):
		RC.errorOut(7, "You have entered a tree type (-t) of 's' but there are labels in your tree that appear more than once!");
	if globs.spec_type == 'm' and any(tips.count(tip) not in [1,2] for tip in tips):
		RC.errorOut(8, "You have entered a tree type (-t) of 'm', species in your tree should appear exactly once or twice.");
	# Some error checking based on the tip labels in the tree.

	if globs.spec_type == 'm':
		hybrid_spec = list(set([tip for tip in tips if tips.count(tip) != 1]));
		for h in hybrid_spec:
			spec_tree = spec_tree.replace(h, h+"*", 1);
	# If the user entered a MUL-tree, some internal re-labeling must be done to those labels that appear twice.

	try:
		sinfo, st = RT.treeParse(spec_tree);
		# Parsing of the species tree.
	except:
		RC.errorOut(9, "Error reading species tree!");
	# Reading the species tree file.

	if globs.label_opt:
		if globs.v != -1:
			print();
			print("# The input species tree with internal nodes labeled:");
			print(st + "\n");
		RC.endProg(starttime);
	# The output if --labeltree is set.

	return sinfo, st;

#############################################################################

def hInParse(sinfo, st, h1_input, h2_input):
	if globs.spec_type == 's':
		hybrid_clades, hybrid_nodes = getHClades(h1_input, sinfo, "h1");
		copy_clades, copy_nodes = getHClades(h2_input, sinfo, "h2");
	# If the input tree is singly-labeled, use the input info from -h1 and -h2 to get the hybrid clades and nodes.

	elif globs.spec_type == 'm':
		mul_copy_clade = [n for n in sinfo if sinfo[n][2] == 'tip' and '*' in n];
		mul_hybrid_clade = [n.replace("*","") for n in mul_copy_clade];

		mul_hybrid_node, mul_hybrid_mono = RT.LCA(mul_hybrid_clade, sinfo);
		mul_copy_node, mul_copy_mono = RT.LCA(mul_copy_clade, sinfo);

		if not mul_hybrid_mono or not mul_copy_mono:
			RC.errorOut(13, "All hybrid clades specified in your MUL-tree must be monophyletic! Hybrid clade identified as: " + ",".join(mul_copy_clade));

		hybrid_clades, hybrid_nodes, copy_clades, copy_nodes = [mul_hybrid_clade], [mul_hybrid_node], [mul_copy_clade], [mul_copy_node];
	# If the input tree is a MUL-tree, we have to determine what the hybrid clades and nodes are.

	return hybrid_clades, hybrid_nodes, copy_clades, copy_nodes;
# Parses the input h nodes.

#############################################################################

def getHClades(h_list, sinfo, h_type):
# This function takes a list of lists of -h1 or -h2 inputs and determines if they are clades or node labels. It then retrieves
# the complete lists of hybrid clades and nodes.
	if h_list:
		if " " in h_list:
			h_clades = h_list.split(" ");
			h_clades = list(map(set, [tmp_h.split(",") for tmp_h in h_clades]));
		else:
			h_clades = list(map(set, [h_list.split(",")]));
		# Split up the input info. If there is a space, multiple nodes/clades have been specified.


		if not all(h in sinfo for hybrid_list in h_clades for h in hybrid_list if not h.isdigit()):
			RC.errorOut(10, "Not all -" + h_type + " species are present in your species tree!");
		if not all("<" + h + ">" in sinfo for hybrid_list in h_clades for h in hybrid_list if h.isdigit()):
			RC.errorOut(11, "Not all -" + h_type + " nodes are present in your species tree!");
		# Some error checking to make sure everything the user input is actually in the tree.

		h_nodes = [];
		for hybrid_clade in h_clades:
			hybrid_clade = list(hybrid_clade);
			if hybrid_clade[0].isdigit():
				h_node = "<" + hybrid_clade[0] + ">";
			# If the input was an internal node, add it to the node list here.
			else:
				h_node, h_mono = RT.LCA(hybrid_clade, sinfo);
				if not h_mono:
					RC.errorOut(12, "All hybrid clades specified h1 and h2 must be monophyletic!");
			# If the input was a clade, retrieve the ancestral node and check if it is monophyletic here.
			if h_node not in h_nodes:
				h_nodes.append(h_node);
			# Add the hybrid node to the nodes list.
	# If the user input anything as -h1 or -h2 this parses it.

	else:
		h_nodes = list(sinfo.keys());
		h_clades = [RT.getClade(node, sinfo) for node in h_nodes];
	# If the user did not specify -h1 or -h2, this adds all possible nodes to the list.

	return h_clades, h_nodes;

#############################################################################

