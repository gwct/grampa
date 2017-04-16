import re, recontree as RT, reconcore as RC, global_vars as globs

#############################################################################

def genMULTrees(hybrid_nodes, copy_nodes, st, sinfo, starttime):
# genMULTree generates all MUL-trees given a species tree and a set of h1 and h2 nodes.
	mul_trees = {};
	mul_num = 1;
	gt_groups = {};

	if globs.lca_opt != 0:
		mul_trees[0] = [st, sinfo, "", "", "", ""];

	if globs.spec_type == 's':
		for hybrid_node in hybrid_nodes:
			hybrid_clade = set(RT.getClade(hybrid_node, sinfo));

			for copy_node in copy_nodes:
				mt_unlabel = buildMultree(hybrid_node, copy_node, st, sinfo);		
				# Building the MUL-tree by passing the species tree to the buildMultree function
				# Input is one node at which to copy the subtree (hybrid_node)
				# and one node at which to place the copy (copy_node)

				if mt_unlabel == "NULL":
					#print("\n*** Warning: H2 node (" + copy_node + ") within hybrid subtree rooted at H1 node (" + hybrid_node + ") or at root of tree. Not building MUL-tree for this combination.\n");
					#print("# ---------------------------");
					continue;
				# Nodes within the hybrid subtree cannot be copy nodes.

				minfo, mt = RT.treeParse(mt_unlabel);
				# Reading the MUL-tree as usual with treeParse.

				if globs.mul_opt:
					outline = hybrid_node + "\t" + copy_node + "\t" + mulPrint(mt, hybrid_clade);
					RC.printWrite(globs.outfilename, globs.main_v, outline);
					continue;

				mul_trees[mul_num] = [mt, minfo, hybrid_clade, hybrid_node, copy_node];
				mul_num += 1;
			
	# This block builds the MUL-trees and prepares the main mul_dict:
	# mul_trees -> mul_num : [MUL-tree string, MUL-tree dict, hybrid clade species, hybrid node, copy node, an empty list to add gene tree groups to later]
	# mul_results -> mul_num : [an initial cumulative score of 0, { gene_num -> [[score for gt/mt combo, # dups for gt/mt combo, # losses for gt/mt combo, maps, node_dups, node_loss],[ditto if multiple maps with same score]]}

	if globs.mul_opt:
		RC.endProg(starttime);
	## If --mulopt is set, this just prints out the MUL-trees and exits.

	return mul_trees;

#############################################################################

def buildMultree(h,p,tree,tree_info):
# This function builds a MUL-tree from a normally labeled species tree.
# Input:
#		h - a node that will have its subtree copied into the MUL-tree
#		p - a node at which the copied subtree will be placed
#		tree - the original species tree WITH internal nodes labeled
#		tree_info - the tree info dictionary (returned by treeparse) from the original tree

	if tree_info[h][2] == 'tip':
		hybrid = h;
	else:
		hybrid = RT.getSubtree(h,tree);
		hybrid = re.sub('<[\d]+>','',hybrid);
	# Gets the subtree of the hybrid node

	if tree_info[p][2] == 'tip':
		copy = p;
	else:
		copy = RT.getSubtree(p,tree);
		copy = re.sub('<[\d]+>','',copy);
	# Gets the subtree of the copy node.

	if (copy in hybrid and copy != hybrid):# or tree_info[p][3] == 'root':
		return "NULL";
	# Copy nodes shouldn't be within the hybrid subtree and shouldn't be at the root.

	if tree_info[h][2] != 'tip':
		for node in tree_info:
			if node in hybrid and tree_info[node][2] == 'tip':
				tree = tree.replace(node, node+"_1");
				copy = copy.replace(node, node+"_1")
	# Some re-labeling if necessary.
	tree = re.sub('<[\d]+>','',tree);

	if tree_info[h][2] == 'tip':
			mul_clade = "(" + copy + "," + hybrid + "*)";
	else:
			mul_clade = "(" + copy + "," + hybrid + ")";

	mul_tree = tree.replace(copy,mul_clade);
	# Combines the clades and replaces the copy clade in the original tree to create the MUL-tree.

	mul_tree = re.sub('<[\d]+>','',mul_tree);
	if tree_info[h][2] != 'tip':
		for node in tree_info:
			if node in hybrid and tree_info[node][2] == 'tip':
				mul_tree = mul_tree.replace(node, node+"*");
		mul_tree = mul_tree.replace("*_1", "");
	# Some relabling of the hybrid species.

	return mul_tree;

#############################################################################

def mulPrint(mul_tree, hybrid_clade):
# For a given MUL-tree, relabels the hybrid clade species to distinguish within tree viewers.

	for spec in hybrid_clade:
		mul_tree = re.sub(spec + '(?!=\*)', spec + '+', mul_tree);
		mul_tree = mul_tree.replace("+*", "*");
	return mul_tree;

#############################################################################