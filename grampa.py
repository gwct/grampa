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
# Multiprocessing implemented Mar/Apr 2017
#############################################################################

import os, timeit, multiprocessing as mp, lib.reconcore as RC, lib.mul_recon as ALG,\
	lib.opt_parse as OP, lib.mul_tree as MT, lib.spec_tree as ST, lib.gene_tree as GT,\
	lib.global_vars as globs, lib.mul_out as OUT
from functools import partial
import cPickle as pickle

def grampa(starttime):

	###########################
	### Prep

	spec_tree_input, gene_tree_input, h1_input, h2_input = OP.optParse(0);
	# Getting the input parameters from optParse.

	if globs.stats:
		import psutil
		pids = [psutil.Process(os.getpid())];		
		prog_start_time = RC.report_stats(outdir=globs.output_directory, stat_start=True);
		step_start_time = prog_start_time;
		globs.v = -1
	# Initializing if --stats is set.

	globs.main_v = -1 if globs.v in [-2,-1] else 1;
	# This is for the tests script and --stats, so nothing prints to the screen.

	if not globs.label_opt:
		RC.printWrite(globs.outfilename, globs.main_v, "# =========================================================================");
		RC.printWrite(globs.outfilename, globs.main_v, "#\t\t  GRAMPA: MUL-tree reconciliation");
		RC.printWrite(globs.outfilename, globs.main_v, "#\t\t\t" + RC.getDateTime());
		if globs.num_procs == 1:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Number of processes:", str(globs.num_procs) + " (serial version)", globs.pad);
		else:
			RC.printWrite(globs.outfilename, globs.main_v, "# LOG " + RC.getLogTime() + ": Number of processes:", str(globs.num_procs), globs.pad);
			pool = mp.Pool(processes = globs.num_procs);
			if globs.stats:
				for result in pool.imap_unordered(RC.getSubPID, range(globs.num_procs)):
					pids.append(result);
		# Prep the pool of processes if the number of processes > 1.
	# Some log info.

	###########################
	### Species tree

	step = RC.printStep(1, "# STEP 1: Reading species tree.");
	sinfo, st = ST.readSpecTree(spec_tree_input);
	hybrid_nodes, copy_nodes = "", "";
	if globs.stats:
		step_start_time = RC.report_stats("Read species trees", pids, step_start_time, prog_start_time, globs.output_directory);

	###########################
	### Hybrid clades

	if globs.lca_opt != 1:
		step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": Parsing hybrid clades.");
		hybrid_clades, hybrid_nodes, copy_clades, copy_nodes = ST.hInParse(sinfo, st, h1_input, h2_input);
		if globs.stats:
			step_start_time = RC.report_stats("Got H nodes", pids, step_start_time, prog_start_time, globs.output_directory);

	###########################
	### Gene trees

	if not globs.mul_opt:
		step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": Reading gene tree(s).");
		gene_trees, gene_trees_filtered, num_skipped = {}, {}, 0;
		# gene_trees_filtered will be a list of lists. One list for each gene tree. If the gene tree passes all
		# filters, it will be [gene tree string, gene tree dict]. Otherwise it will be [Filter message].

		if os.path.isfile(gene_tree_input):
			try:
				gene_num = 1;
				for line in open(gene_tree_input, "r"):
					gene_trees[gene_num] = line.strip();
					gene_num += 1;
			except:
				RC.errorOut(13, "Error reading gene trees file!");
		# If the gene tree input is a file, try reading the file.
		else:
			gene_trees[1] = gene_tree_input;
		# Otherwise try it as a tree string.

		if globs.num_procs != 1:
			for result in pool.imap_unordered(GT.readGeneTree, gene_trees.iteritems()):
				gene_trees_filtered[result[0]] = result[1];
				if result[2]:
					num_skipped += 1;	
		### PARALLEL VERSION

		else:
			for gene_num, gene_tree in gene_trees.iteritems():
				gn, gt, skip_flag = GT.readGeneTree((gene_num, gene_tree));
				gene_trees_filtered[gene_num] = gt;
				if skip_flag:
					num_skipped += 1;
		### SERIAL VERSION
		# Parsing the gene trees to get info about each node in each tree.

		if num_skipped == len(gene_trees):
			RC.errorOut(14, "Couldn't find any gene trees in your gene tree input file (-g)!");
		# If the input is a file, we assume each line contains one gene tree.
	else:
		gene_trees_filtered = "";
	## Reading the input files.
	if globs.stats:
		step_start_time = RC.report_stats("Read gene trees", pids, step_start_time, prog_start_time, globs.output_directory);

	###########################

	OUT.logOut(st, str(len(gene_trees_filtered)), hybrid_nodes, copy_nodes, gene_tree_input, h1_input, h2_input);
	# Printing out a bunch of log info.

	###########################
	### Build MUL-trees

	if globs.lca_opt == 1 and globs.spec_type == 's':
	# If the user specifies to do reconciliations ONLY to the singly-labeled tree, then the mul_trees dict should
	# have only that tree in it. Also, the collapseGroups step can be skipped.
		mul_trees = { 0 : [st, sinfo, "", "", ""] };
		gene_trees = gene_trees_filtered;

	elif globs.spec_type == 's':
		step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": Building " + str(len(hybrid_nodes) * len(copy_nodes)) + " MUL-tree(s).");
		mul_trees = MT.genMULTrees(hybrid_nodes, copy_nodes, st, sinfo, starttime);
		# genMULTree does all the parsing of the h1 and h2 clades and calls the buildMulTree function.

	elif globs.spec_type == 'm':
		mul_trees = { 1 : [st, sinfo, hybrid_clades[0], hybrid_nodes[0], copy_nodes[0], []] };
	# If the input tree is a MUL-tree, we just need to set some variables.

		if globs.stats:
			step_start_time = RC.report_stats("Built MUL-trees", pids, step_start_time, prog_start_time, globs.output_directory);
	# The two main MUL-tree variables are described here:
	# mul_trees -> mul_num : [MUL-tree string, MUL-tree dict, hybrid clade species, hybrid node, copy node, an empty list to add gene tree groups to later]

	###########################
	### Collapse groups
	if globs.lca_opt != 1:
		step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": Checking gene tree groupings.");

		if globs.num_procs != 1 and globs.spec_type == 's':
			for result in pool.imap_unordered(partial(ALG.collapseGroups, gene_trees_filtered_cg=gene_trees_filtered, spec_type_cg=globs.spec_type, pickle_dir=globs.pickle_dir), mul_trees.iteritems()):
				pass;
		### PARALLEL AND SINGLY LABELED TREE VERSION

		else:
			for mul_num, mul_tree in mul_trees.iteritems():
				ALG.collapseGroups((mul_num, mul_tree), gene_trees_filtered, globs.spec_type, globs.pickle_dir);
		### SERIAL AND MUL-TREE VERSION
		# If the input tree is singly-labeled, do some multi-processing, else if it is a MUL-tree, just run the function.

	num_skipped = OUT.checkOut(mul_trees, num_skipped, gene_trees_filtered);
	gene_trees, step = OUT.filterOut(num_skipped, step, gene_trees_filtered)
	
	# Write the filtered trees to a file, or not if no filtering was done.
	if globs.check_nums:
		RC.endProg(starttime);
	# If --checknums is set, exit the program here.
	RC.printWrite(globs.outfilename, globs.main_v, "# ----------------------------------------");

	if globs.stats:
		step_start_time = RC.report_stats("Collapsed groups", pids, step_start_time, prog_start_time, globs.output_directory);

	###########################
	### Reconciliations

	step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": Beginning reconciliations.");
	if globs.stats:
		step_start_time = RC.report_stats("Pre-recon", pids, step_start_time, prog_start_time, globs.output_directory);
		mt_start_time = step_start_time;
	
	all_scores, min_num, min_score, min_maps = {}, '', 9999999, {};

	if globs.num_procs != 1:
		for result in pool.imap_unordered(partial(ALG.mulRecon, gene_trees=gene_trees, v=globs.v, pickle_dir=globs.pickle_dir), mul_trees.iteritems()):
			all_scores[result[0]] = result[1];
			if result[1] < min_score:
				min_num = result[0];
				min_score = result[1];
			if globs.stats:
				mt_start_time = RC.report_stats("MT-" + str(result[0]), pids, mt_start_time, prog_start_time, globs.pickle_dir);
			del result;
		pool.terminate();
		if globs.stats:
			pids = [pids[0]];
	### PARALLEL VERSION

	else:
		for mul_num, mul_tree in mul_trees.iteritems():
			mn, cur_score = ALG.mulRecon((mul_num, mul_tree), gene_trees, globs.v, globs.pickle_dir);
			all_scores[mul_num] = cur_score;
			if cur_score < min_score:
				min_num = mul_num;
				min_score = cur_score;
			if globs.stats:
				mt_start_time = RC.report_stats("MT-" + str(mul_num), pids, step_start_time, prog_start_time, globs.output_directory);
	### SERIAL VERSION

	if globs.stats:
		step_start_time = RC.report_stats("Post-recon", pids, step_start_time, prog_start_time, globs.output_directory);

	###########################
	### Output

	min_tree = mul_trees[min_num];
	min_maps = ALG.mulRecon((min_num, min_tree), gene_trees, globs.v, globs.pickle_dir, retmap=True);
	## This gets the mul_tree with the min score along with its maps.

	multiple_maps = OUT.detOut(gene_trees, min_tree, min_num, min_maps);
	OUT.mainOut(mul_trees, all_scores, min_num, min_score, min_maps, multiple_maps);

	if globs.stats:
		step_start_time = RC.report_stats("Output done", pids, step_start_time, prog_start_time, globs.output_directory);
		print("# " + "-" * 100);

	###########################
	### Begin orthology prediction block. *BETA*

	if globs.orth_opt:
		import lib.orth_label as OL
		step = RC.printStep(step, "# STEP " + str(step) + " " + RC.getLogTime() +  ": *BETA* Mapping polyploid genes.");
		min_tree = mul_trees[min_num][0];
		min_clade = mul_trees[min_num][2];
		OL.orthLabel(gene_trees, min_maps, min_tree, min_clade);

	###########################

	return;

#############################################################################

if __name__ == '__main__':
# Necessary for multiprocessing to work on Windows.
	start = timeit.default_timer();
	globs.init();
	grampa(start);
	print("# LOG " + RC.getLogTime() + ": Done!");
	RC.endProg(start);

#############################################################################







	# if globs.v == 0:
	# 	numiters = len(mul_dict) * len(gene_trees_filtered);
	# 	numbars = 0;
	# 	donepercent = [];
	# 	itercount = 0;
	# Stuff for the loading bar...

	# Check if the current MUL-tree has a lower score than the previous lower score. If so, save it as the lowest_score.
	# if globs.v == 0 and numiters > 100:
	# 	pstring = "100.0% complete\n";
	# 	sys.stderr.write('\b' * len(pstring) + pstring);