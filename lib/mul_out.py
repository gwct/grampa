import os, lib.reconcore as RC, lib.mul_tree as MT, lib.global_vars as globs
import pickle

#############################################################################

def logOut(st, gene_trees_filtered, hybrid_nodes, copy_nodes, gene_tree_input,  h1_input, h2_input):
	RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": The input species tree with internal nodes labeled:", st, globs.pad);
	if globs.spec_type == 's' and globs.lca_opt != 1:
		RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Input species tree is:", "Singly-labeled", globs.pad);
		if not h1_input:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": No H1 node defined", "Searching all possible H1 nodes.", globs.pad);
		else:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": H1 node(s) identified as:", ",".join(hybrid_nodes), globs.pad);
		if not h2_input:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": No H2 node defined", "Searching all possible H2 nodes.", globs.pad);
		else:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": H2 node(s) identified as:", ",".join(copy_nodes), globs.pad);
	elif globs.spec_type == 'm':
		RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Input species tree is:", "MUL-tree", globs.pad);
	if not globs.mul_opt:
		if os.path.isfile(gene_tree_input):
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Using gene trees in file:", gene_tree_input, globs.pad);
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Total gene trees:", gene_trees_filtered, globs.pad);
		else:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": The input gene tree with internal nodes labeled:", gene_tree_input, globs.pad);
	RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Main results and logs will be written to file:", globs.outfilename, globs.pad);
	if globs.lca_opt == 1:
		RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Reconciling to:", "Singly-labeled tree only.", globs.pad);
	else:
		if not globs.mul_opt:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": The number of groups for each tree will be calculated:", globs.checkfilename, globs.pad);
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Filtered trees will be saved (if necessary):", globs.gene_file_filtered, globs.pad);
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Filtering gene trees with # of polyploid groups over:", str(globs.cap), globs.pad);
			if not globs.check_nums:
				RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Detailed results will be written to file:", globs.detoutfilename, globs.pad);
				if globs.maps_opt:
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Detailed output will contain:", "Reconciliation scores and maps", globs.pad);
				else:
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Detailed output will contain:", "Reconciliation scores only", globs.pad);
				if globs.lca_opt == 0:
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Reconciling to:", "MUL-trees only", globs.pad);
				elif globs.lca_opt == 2:
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Reconciling to:", "Singly-labeled and MUL-trees", globs.pad);
				if globs.orth_opt:
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": *BETA* Attempting to discern relationships of polyploid genes (paralog vs. homoeolog).")
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": *BETA* Writing relationships to:", globs.orth_file_name, globs.pad);
					RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": *BETA* Writing labeled tres to:", globs.labeled_tree_file, globs.pad);
			elif globs.check_nums:
				RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": --checknums set. NOT doing reconciliations, just running some numbers for you.");
		elif globs.mul_opt:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": --multree set. NOT doing reconciliations, just building your MUL-trees.");

#############################################################################

def checkOut(mul_trees, num_skipped, gene_trees_filtered):
	checkfile = open(globs.checkfilename, "a");
	for mul_num in mul_trees:
		if mul_num == 0:
			continue;
		if globs.spec_type == 's':
			checkfile.write("# MT-" + str(mul_num) + ":" + MT.mulPrint(mul_trees[mul_num][0], mul_trees[mul_num][2]) + "\tH1 Node:" + mul_trees[mul_num][3] + "\tH2 Node:" + mul_trees[mul_num][4] + "\n");
		elif globs.spec_type == 'm':
			checkfile.write("# MT-" + str(mul_num) + ":" + MT.mulPrint(mul_trees[mul_num][0], mul_trees[mul_num][2]) + "\n");

		groupfilename = os.path.join(globs.pickle_dir, str(mul_num) + "_groups.pickle");
		cur_groups = pickle.load(open(groupfilename, "rb"));

		for gene_num in gene_trees_filtered:
			if len(gene_trees_filtered[gene_num]) == 1:
				continue;
			gt_groups, gt_fixed = cur_groups[gene_num][0], cur_groups[gene_num][1];
			outline = "GT-" + str(gene_num) + " to MT-" + str(mul_num) + "\t";
			num_groups = len(gt_groups);
			num_fixed = len(gt_fixed);
			outline += str(num_groups) + "\t" + str(num_fixed) + "\t" + str(2**num_groups);
			if num_groups > globs.cap:
				gene_trees_filtered[gene_num] = ["# Number of groups over group cap (-c set to " + str(globs.cap) + ") -- Filtering."];
				outline += "\tNumber of groups over group cap (-c set to " + str(globs.cap) + ") -- Filtering.";
				num_skipped += 1;
			checkfile.write(outline + "\n");
		checkfile.write("# ----------------------------------\n");
	checkfile.close();

	return num_skipped;
	# This block handles output to the checknums file and filters any trees over the group cap (-c).
	# The call of the important collapseGroups function that groups polyploid clades in the gene trees to speed up the reconciliations to MUL-trees.

#############################################################################

def filterOut(num_skipped, step, gene_trees_filtered):
	gene_trees = {};
	if num_skipped != 0:
		RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Filtered trees:", str(num_skipped), globs.pad);
		RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Final tree count for reconciliations:", str(len(gene_trees_filtered) - num_skipped), globs.pad);
		step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": Writing filtered gene trees to file.");
		filtered_file = open(globs.gene_file_filtered, "w");
		for gene_num in gene_trees_filtered:
			filtered_file.write(gene_trees_filtered[gene_num][0] + "\n");
			if len(gene_trees_filtered[gene_num]) != 1:
				gene_trees[gene_num] = gene_trees_filtered[gene_num]
		filtered_file.close();
	else:
		RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": No trees filtered! Using your original set.");
		gene_trees = gene_trees_filtered;

	return gene_trees, step;

#############################################################################

def mainOut(mul_trees, all_scores, min_num, min_score, min_maps, multiple_maps):
	RC.printWrite(globs.outfilename, globs.v, "# Tree #\tH1 node\tH2 node\tTree string\tTotal score");
	for mul_num, mul_tree in mul_trees.iteritems():
		if mul_num == min_num:
			min_tree = mul_tree;

		mt, hybrid_clade, hybrid_node, copy_node = mul_tree[0], mul_tree[2], mul_tree[3], mul_tree[4];
		cur_score = all_scores[mul_num];
		if mul_num == 0:
			outline = "ST\t\t" + mt + "\t" + str(cur_score);
		elif globs.spec_type == 's':
			outline = "MT-" + str(mul_num) + "\t" + hybrid_node + "\t" + copy_node + "\t" + MT.mulPrint(mt, hybrid_clade) + "\t" + str(cur_score);
		elif globs.spec_type == 'm':
			outline = "MT-" + str(mul_num) + "\t" + MT.mulPrint(mt, hybrid_clade) + "\t" + str(cur_score);
	# Save the main output for ALL MUL-trees.
		RC.printWrite(globs.outfilename, globs.v, outline);
	RC.printWrite(globs.outfilename, globs.main_v, "# ----------------------------------------");
	if globs.spec_type == 's':
		if min_num != 0:
			RC.printWrite(globs.outfilename, globs.v, "The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" \
				+ MT.mulPrint(min_tree[0], min_tree[2]));
		else:
			RC.printWrite(globs.outfilename, globs.v, "The tree with the minimum parsimony score is the singly-labled tree (ST):\t" + min_tree[0]);
		RC.printWrite(globs.outfilename, globs.v, "Score = " + str(min_score));
		RC.printWrite(globs.outfilename, globs.main_v, "# ----------------------------------------");
		RC.printWrite(globs.outfilename, globs.v, "Species tree node\t# dups mapped");

		mt, min_dict = min_tree[0], mul_tree[1];
		from collections import defaultdict
		main_dups = defaultdict(int);
		for m in min_maps:
			maps, dups = min_maps[m][0][3], min_maps[m][0][4];
			for gt_node in dups:
				if dups[gt_node] != 0:
					main_dups[maps[gt_node][0]] += 1;
		for node in main_dups:
			RC.printWrite(globs.outfilename, globs.v, node + "\t" + str(main_dups[node]));
	if multiple_maps != 0:
		num_str = " trees have ";
		if multiple_maps == 1:
			num_str = " tree has ";
		multiple_outline = "# MSG: " + str(multiple_maps) + " gene" + num_str + "multiple maps to the species tree with equal scores. Only one of these maps is (randomly) chosen in the final duplication counts. See detailed output file for more info."
		RC.printWrite(globs.outfilename, globs.main_v, "# ----------------------------------------");
		RC.printWrite(globs.outfilename, globs.main_v, multiple_outline)
	# Output to the main file.

#############################################################################

def detOut(gene_trees, min_tree, min_num, min_maps):
	import gene_tree as GT
	multiple_maps = 0;
	det_header = "# GT/MT combo\t# dups\t# losses\tTotal score";
	if globs.maps_opt:
		det_header += "\tMaps"
	RC.printWrite(globs.detoutfilename, globs.v, det_header);
	# Header info for the detailed output file.

	for gene_num, cur_maps in min_maps.iteritems():
		mul_map_string = "";
		if len(cur_maps) != 1:
			multiple_maps += 1;
			RC.printWrite(globs.detoutfilename, globs.v, "* GT-" + str(gene_num) + " to MT-" + str(min_num) + "\t" + str(len(cur_maps)) + " maps found!");
			mul_map_string = "* "

		for cur_map in cur_maps:
			outline = mul_map_string + "GT-" + str(gene_num) + " to MT-" + str(min_num) + "\t";
			outline += str(cur_map[1]) + "\t" + str(cur_map[2]) + "\t" + str(cur_map[0]);
			
			if globs.maps_opt:
				outline += "\t" + GT.detailedOut(gene_trees[gene_num][0], gene_trees[gene_num][1], cur_map[3], cur_map[4], cur_map[5]);
			RC.printWrite(globs.detoutfilename, globs.v, outline);

	return multiple_maps;
	# Output to the detailed file.

#############################################################################