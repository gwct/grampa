import grampa_lib.recontree as RT
import grampa_lib.reconcore as RC

#############################################################################
def readGeneTree(gene_tree_input):
# Reading gene tree dictionaries.
	gene_num, gene_tree = gene_tree_input;

	if gene_tree.strip() == '':
		gene_tree = ["# Empty line -- Filtering."];
		return gene_num, gene_tree, True;
	# Handles empty lines in the gene tree file.

	try:
		gene_tree = RT.remBranchLength(gene_tree);
		ginfo, gt = RT.treeParse(gene_tree);
		gene_tree = [gt,ginfo];
	except:
		gene_tree = ["# Error reading this tree! -- Filtering."];
		return gene_num, gene_tree, True;
	# Tries the gene tree parsing code and if anything goes wrong, catches exception and filters the tree.

	num_tips = len([n for n in ginfo if ginfo[n][2] == 'tip']);
	num_internal = len([n for n in ginfo if ginfo[n][2] != 'tip']);

	if num_internal != num_tips - 1:
		gene_tree = ["# This line may not contain a tree, or if so it may be unrooted or contain a polytomy -- Filtering."];
		return gene_num, gene_tree, True;
	# Another check for gene tree parsing and formatting errors.

	return gene_num, gene_tree, False;

#############################################################################

def detailedOut(gt, ginfo, maps, dups, losses):
# Outputs the gene trees for the --maps option. Adds the maps, dups, and losses to each node label.
	import re
	for node in ginfo:
		cur_map = maps[node][0];
		if "*" not in cur_map:
			cur_map += "+";
		node_string = node + "[" + cur_map + "-" + str(dups[node]) + "]";
		gt = re.sub("(?<![\[])" + node, node_string, gt);	

	return gt;

#############################################################################