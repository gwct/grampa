import os, recontree as RT, reconcore as RC, optparse as OP

#############################################################################

def readGeneTree(gene_tree_input):
# Reading the gene trees file.
	gene_trees_filtered = [];
	# gene_trees_filtered will be a list of lists. One list for each gene tree. If the gene tree passes all
	# filters, it will be [gene tree string, gene tree dict]. Otherwise it will be [Filter message].

	num_skipped = 0;
	if os.path.isfile(gene_tree_input):
		try:
			gene_trees = open(gene_tree_input, "r").readlines();
		except:
			RC.errorOut(13, "Error reading gene trees file!");
		# Try reading the file.

		gene_num = -1;

		for gene_tree in gene_trees:
			gene_num += 1;
			if gene_tree.strip() == '':
				gene_trees_filtered.append(["# Empty line -- Filtering."]);
				num_skipped += 1;
				continue;
			# Handles empty lines in the gene tree file.

			try:
				gene_tree = RT.remBranchLength(gene_tree);
				ginfo, gt = RT.treeParse(gene_tree);
				gene_trees_filtered.append([gt,ginfo]);
			except:
				gene_trees_filtered.append(["# Error reading this tree! -- Filtering."]);
				num_skipped += 1;
				continue;
			# Tries the gene tree parsing code and if anything goes wrong, catches exception and filters the tree.

			if len([g for g in ginfo if ginfo[g][2] != 'tip']) != len([g for g in ginfo if ginfo[g][2] == 'tip']) - 1:
				gene_trees_filtered[gene_num] = ["# This line may not contain a tree, or if so it may be unrooted -- Filtering."]
				num_skipped += 1;
				continue;
			# Another check for gene tree parsing and formatting errors.

		if num_skipped == len(gene_trees):
			RC.errorOut(14, "Couldn't find any gene trees in your gene tree input file (-g)!");
	# If the input is a file, we assume each line contains one gene tree.

	else:
		gene_tree = gene_tree_input;
		gt_err_flag = 0;

		gene_tree = RT.remBranchLength(gene_tree);
		ginfo, gt = RT.treeParse(gene_tree);
		gene_trees_filtered.append([gt,ginfo]);

		if len([g for g in ginfo if ginfo[g][2] != 'tip']) != len([g for g in ginfo if ginfo[g][2] == 'tip']) - 1:
			gt_err_flag = 1;

		if gt_err_flag == 1:
			RC.errorOut(15, "Error reading your input gene tree string! Make sure it's rooted!")
	# If the input string is a filename, read the file. Otherwise, just try it as a newick string.

	return gene_trees_filtered, num_skipped;

#############################################################################

def detailedOut(maps, dups, losses, v, detoutfilename):
	for node in maps:
		outline = "\t" + node + "\t" + maps[node][0] + "\t" + str(dups[node]) + "\t" + str(losses[node]);
		RC.printWrite(detoutfilename, v, outline);

#############################################################################