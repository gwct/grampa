#!/usr/bin/python
#############################################################################
# Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.
# This is the main interface and handles user options and tree searching.
#
# Gregg Thomas
# Fall 2015
# Combo algorithm implemented Spring 2016
# Filter update implemented Feb. 2017
# Orthology labeling started in Mar. 2017 (BETA)
#############################################################################

import sys, os, re, time, lib.recontree as RT, lib.reconcore as RC, lib.mul_recon as ALG, lib.lca_check as LCHECK, \
lib.opt_parse as OP, lib.mul_tree as MT, lib.std_lca as STD, lib.spec_tree as ST, lib.gene_tree as GT

#############################################################################
starttime = time.time();

spec_tree_input, spec_type, gene_tree_input, h1_input, h2_input, cap, outfilename, \
checkfilename, detoutfilename, gene_file_filtered, labeled_tree_file, orth_file, v, lca_opt, label_opt, \
check_nums, mul_opt, orth_opt, maps_opt = OP.optParse(0);
# Getting the input parameters from optParse.

main_v = 1 if v!= -1 else -1;
if not label_opt:
	RC.printWrite(outfilename, main_v, "# =========================================================================");
	RC.printWrite(outfilename, main_v, "#\t\t\tGRAMPA: MUL-tree reconciliation");
	RC.printWrite(outfilename, main_v, "#\t\t\t" + RC.getDateTime());

step = RC.printStep(1, "# STEP 1: Reading species tree.", v);
sinfo, st = ST.readSpecTree(spec_tree_input, spec_type, label_opt, v);
step = RC.printStep(step, "# STEP " + str(step) + ": Parsing hybrid clades.", v);
hybrid_clades, hybrid_nodes, copy_clades, copy_nodes = ST.hInParse(sinfo, st, h1_input, h2_input, spec_type);
if not mul_opt:
	step = RC.printStep(step, "# STEP " + str(step) + ": Reading gene tree(s).", v);
	gene_trees_filtered, num_skipped = GT.readGeneTree(gene_tree_input);

lowest_score = ["ST", st, sinfo, "", "", "", 9999999];
# Lowest score will hold the main outputs from the lowest scoring tree: tree type, tree, tree info, hybrid_clade, maps for each gene tree, duplication nodes, score;

#############################################################################
### Begin input info block!
pad = 65
RC.printWrite(outfilename, main_v, "# LOG: The input species tree with internal nodes labeled:", st, pad);
if spec_type == 's' and lca_opt != 1:
	RC.printWrite(outfilename, main_v, "# LOG: Input species tree is:", "Singly-labeled", pad);
	if not h1_input:
		RC.printWrite(outfilename, main_v, "# LOG: No H1 node defined", "Searching all possible H1 nodes.", pad);
	else:
		RC.printWrite(outfilename, main_v, "# LOG: H1 node(s) identified as:", ",".join(hybrid_nodes), pad);
	if not h2_input:
		RC.printWrite(outfilename, main_v, "# LOG: No H2 node defined", "Searching all possible H2 nodes.", pad);
	else:
		RC.printWrite(outfilename, main_v, "# LOG: H2 node(s) identified as:", ",".join(copy_nodes), pad);
elif spec_type == 'm':
	RC.printWrite(outfilename, main_v, "# LOG: Input species tree is:", "MUL-tree", pad);
if not mul_opt:
	if os.path.isfile(gene_tree_input):
		RC.printWrite(outfilename, main_v, "# LOG: Using gene trees in file:", gene_tree_input, pad);
		if lca_opt != 1:
			RC.printWrite(outfilename, main_v, "# LOG: Filtering gene trees with # of polyploid groups over:", str(cap), pad);
	else:
		RC.printWrite(outfilename, main_v, "# LOG: The input gene tree with internal nodes labeled:", gene_tree_input, pad);
RC.printWrite(outfilename, main_v, "# LOG: Main results will be written to file:", outfilename, pad);
if lca_opt == 1:
	RC.printWrite(outfilename, main_v, "# LOG: Reconciling to:", "Singly-labeled tree only.", pad);
else:
	if not mul_opt and lca_opt != 1:
		RC.printWrite(outfilename, main_v, "# LOG: The number of groups for each tree will be calculated:", checkfilename, pad);
		RC.printWrite(outfilename, main_v, "# LOG: Filtered trees will be saved (if necessary):", gene_file_filtered, pad);
		if not check_nums:
			RC.printWrite(outfilename, main_v, "# LOG: Detailed results will be written to file:", detoutfilename, pad);
			if maps_opt:
				RC.printWrite(outfilename, main_v, "# LOG: Detailed output will contain reconciliation scores and maps.");
			else:
				RC.printWrite(outfilename, main_v, "# LOG: Detailed output will contain reconciliation scores only.");
			if lca_opt == 0:
				RC.printWrite(outfilename, main_v, "# LOG: Reconciling to:", "MUL-trees only", pad);
			elif lca_opt == 2:
				RC.printWrite(outfilename, main_v, "# LOG: Reconciling to:", "Singly-labeled and MUL-trees", pad);

			if orth_opt:
				RC.printWrite(outfilename, main_v, "# LOG: *BETA* Attempting to discern relationships of polyploid genes (paralog vs. homoeolog).")
				RC.printWrite(outfilename, main_v, "# LOG: *BETA* Writing relationships to:", orth_file, pad);
				RC.printWrite(outfilename, main_v, "# LOG: *BETA* Writing labeled tres to:", labeled_tree_file, pad);
		elif check_nums:
			RC.printWrite(outfilename, main_v, "# LOG: --checknums set. NOT doing reconciliations, just running some numbers for you...");
	elif mul_opt:
		RC.printWrite(outfilename, main_v, "# LOG: --multree set. NOT doing reconciliations, just building your MUL-trees...");
RC.printWrite(outfilename, main_v, "# ---------");
### End input info block!
#############################################################################	
if lca_opt != 1:
	if spec_type == 's':
		step = RC.printStep(step, "# STEP " + str(step) + ": Building " + str(len(hybrid_nodes) * len(copy_nodes)) + " MUL-tree(s).", v);
		mul_dict = MT.genMULTree(hybrid_nodes, copy_nodes, st, sinfo, spec_type, outfilename, mul_opt, lca_opt, v, main_v, starttime);
		# genMULTree does all the parsing of the h1 and h2 clades and calls the buildMulTree function.
	elif spec_type == 'm':
		mul_dict = { 1 : [st, sinfo, hybrid_clades[0], hybrid_nodes[0], copy_nodes[0], 0, []] };
	# mul_dict -> mul_num : [MUL-tree string, MUL-tree dict, hybrid clade species, hybrid node, copy node, an initial score of 0, an empty list to add gene tree groups to later]
## Generate the MUL trees if not doing reconcilition to ST only (ie if lca_opt != 1).
#############################################################################
if lca_opt != 1:
# Do not need to check groups if only doing standard LCA.
	step = RC.printStep(step, "# STEP " + str(step) + ": Checking gene tree groupings.", v);

	mul_dict, gene_trees_filtered, num_skipped = ALG.collapseGroups(mul_dict, gene_trees_filtered, checkfilename, spec_type, num_skipped, cap, v);
	# The call of the important collapseGroups function that groups polyploid clades in the gene trees to speed up the reconciliations to MUL-trees.

	if num_skipped != 0:
		RC.printWrite(outfilename, v, "# LOG: Filtered " + str(num_skipped) + " trees.");
		step = RC.printStep(step, "# STEP " + str(step) + ": Writing filtered gene trees to file.", v);
		filtered_file = open(gene_file_filtered, "w");
		for gt in gene_trees_filtered:
			filtered_file.write(gt[0] + "\n");
		filtered_file.close();
	else:
		RC.printWrite(outfilename, v, "# LOG: No trees filtered! Using your original set.");
	# Write the filtered trees to a file, or not if no filtering was done.

	if check_nums:
		RC.endProg(starttime, outfilename, main_v);
	# If --checknums is set, exit the program here.

#############################################################################
if spec_type == 's' and lca_opt != 0:
	step = RC.printStep(step, "# STEP " + str(step) + ": Doing standard reconciliation to your singly-labeled tree.", v);
	lowest_score = STD.stdLCA(st, sinfo, gene_trees_filtered, outfilename, detoutfilename, lca_opt, maps_opt, v);
	if lca_opt == 1:
		RC.endProg(starttime, outfilename, main_v);

#############################################################################
step = RC.printStep(step, "# STEP " + str(step) + ": Beginning reconciliations to MUL-trees.", v);
if v == 0:
	numiters = len(mul_dict) * len(gene_trees_filtered);
	numbars = 0;
	donepercent = [];
	itercount = 0;
# Stuff for the loading bar...

for mul_num in mul_dict:
	mt, minfo, hybrid_clade, hybrid_node, copy_node, group_list = mul_dict[mul_num][0], mul_dict[mul_num][1], \
	mul_dict[mul_num][2], mul_dict[mul_num][3], mul_dict[mul_num][4], mul_dict[mul_num][6], 

	gt_maps, gt_dups, multiple_maps = [], [], 0;

	tot_node_counts = {};
	for node in minfo:
		tot_node_counts[node] = [0,0];

	if spec_type == 's':
		RC.printWrite(detoutfilename, 0, "MT-" + str(mul_num) + ":" + MT.mulPrint(mt, hybrid_clade) + "\tH1 Node:" + hybrid_node + "\tH2 Node:" + copy_node);
	elif spec_type == 'm':
		RC.printWrite(detoutfilename, 0, "MT-" + str(mul_num) + ":" + MT.mulPrint(mt, hybrid_clade));
	gene_num = -1;

	for gene_tree in gene_trees_filtered:
		if v == 0 and numiters > 100:
			numbars, donepercent = RC.loadingBar(itercount, numiters, donepercent, numbars);
			itercount = itercount + 1;
		# Only the loading bar displays when the program is running if -v is set to 0.

		cur_maps, cur_dups = [], [];

		gene_num += 1;
		if len(gene_tree) == 1:
			gt_maps.append(cur_maps);
			gt_dups.append(cur_dups);
			continue;
		# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

		gt_groups, gt_fixed = group_list[gene_num][0], group_list[gene_num][1];
		gt, ginfo = gene_tree;
		# Retrieve gene tree info and collapsed groups for this gene tree-MUL-tree combo

		gt_results = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, gt_groups, gt_fixed, v);
		# The call of the reconciliation algorithm! On the current gene tree with the current MUL-tree.
		# gt_results: [[score, dup_score, loss_score, maps, node_dups, node_loss],[ditto]]

		mul_dict[mul_num][5] += gt_results[0][0];
		mul_map_string = "";
		if len(gt_results) != 1:
			multiple_maps += 1;
			RC.printWrite(detoutfilename, v, "* GT-" + str(gene_num+1) + " to MT-" + str(mul_num) + "\t" + str(len(gt_results)) + " maps found!");
			mul_map_string = "* "

		for gt_result in gt_results:
			outline = mul_map_string + "GT-" + str(gene_num+1) + " to MT-" + str(mul_num) + "\t";
			outline += str(gt_result[1]) + "\t" + str(gt_result[2]) + "\t" + str(gt_result[0]);
			RC.printWrite(detoutfilename, v, outline);

			if maps_opt:
				GT.detailedOut(gt_result[3], gt_result[4], gt_result[5], v, detoutfilename);

			cur_maps.append(gt_result[3]);
			cur_dups.append(gt_result[4]);

		gt_maps.append(cur_maps);
		gt_dups.append(cur_dups);

	RC.printWrite(detoutfilename, v, "Gene trees with multiple maps:\t" + str(multiple_maps));
	RC.printWrite(detoutfilename, v, "Total parsimony score for MT-" + str(mul_num) + ": " + str(mul_dict[mul_num][5]));
	RC.printWrite(detoutfilename, v, "# ---------------------------");
	if spec_type == 's':
		RC.printWrite(outfilename, 0, "MT-" + str(mul_num) + "\t" + hybrid_node + "\t" + copy_node + "\t" + MT.mulPrint(mt, hybrid_clade) + "\t" + str(mul_dict[mul_num][5]));
	elif spec_type == 'm':
		RC.printWrite(outfilename, 0, "MT-" + str(mul_num) + "\t" + MT.mulPrint(mt, hybrid_clade) + "\t" + str(mul_dict[mul_num][5]));
	# Print the total score for the current MUL-tree.

	if mul_dict[mul_num][5] < lowest_score[6]:
		lowest_score = [mul_num, mt, minfo, hybrid_clade, gt_maps, gt_dups, mul_dict[mul_num][5]];
	# Check if the current MUL-tree has a lower score than the previous lower score. If so, save it as the lowest_score.
	
if v == 0 and numiters > 100:
	pstring = "100.0% complete\n";
	sys.stderr.write('\b' * len(pstring) + pstring);
#############################################################################
score, min_num, min_tree, min_clade = lowest_score[6], lowest_score[0], lowest_score[1], lowest_score[3];
if spec_type == 's':
	if min_num != "ST":
		RC.printWrite(outfilename, 0, "# The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + MT.mulPrint(min_tree, min_clade));
	else:
		RC.printWrite(outfilename, 0, "# The tree with the minimum parsimony score is the singly-labled tree (ST):\t" + st);
	RC.printWrite(outfilename, 0, "# Score = " + str(score));

if orth_opt:
	step = RC.printStep(step, "# STEP " + str(step) + ": *BETA* Mapping polyploid genes...", v);

	orthfile = open(orth_file, "a");
	lfile = open(labeled_tree_file, "a");

	min_maps, min_dups = lowest_score[4], lowest_score[5];

	gene_num = -1;

	orthfile.write(lowest_score[1] + "\n");
	orthfile.write("------------------\n");

	for gene_tree in gene_trees_filtered:
		gene_num += 1;
		if len(gene_tree) == 1:
			continue;
		# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

		cur_maps, cur_dups = min_maps[gene_num], min_dups[gene_num];
		gt, ginfo = gene_tree;
		# Retrieve gene tree info and collapsed groups for this gene tree-MUL-tree combo

		if len(cur_maps) != 1:
			orthfile.write(str(gene_num+1) + "\t* " + str(len(cur_maps)) + " maps tied for lowest score. Mapping orthologies to each set of maps!\n");
			lfile.write(str(gene_num+1) + "\t* " + str(len(cur_maps)) + " maps tied for lowest score. Labeling trees for each set of maps!\n");

		for x in range(len(cur_maps)):
			cur_map, cur_dup = cur_maps[x], cur_dups[x];

			LCHECK.lcaCheck(gt, ginfo, cur_map, cur_dup, min_clade, gene_num, lfile, orthfile);

	lfile.close();
	orthfile.close();
### NEW SECTION TO LABEL HOMEOLOGS/PARALOGS

if v != -1:
	print("# LOG: Done!");
# Final output block.
RC.endProg(starttime, outfilename, main_v);
### End scoring block and end program!

