#!/usr/bin/python
#############################################################################
# Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.
# This is the main interface and handles user options and tree searching.
#
# Gregg Thomas
# Fall 2015, Combo algorithm implemented Spring 2016
#############################################################################

import sys, os, re, time, argparse, lib.recontree as RT, lib.reconcore as RC, lib.mul_recon as ALG

############################################
#Function Definitions
############################################
def optParse(errorflag):
# This function handles the command line options.

	parser = argparse.ArgumentParser(description="GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.");

	parser.add_argument("-s", dest="spec_tree", help="A file containing a bifurcating species tree in newick format on which to search for polyploid events.");
	parser.add_argument("-t", dest="spec_tree_type", help="[m or s] -- m: input species tree is a MUL-tree. s: input species tree is a standard tree. Default: s", default="s");
	parser.add_argument("-g", dest="gene_input", help="A file containing one or more newick formatted gene trees to reconcile. The labels in the gene tree must end with '_[species name]' and contain no other underscores.");
	parser.add_argument("-h1", dest="h1_spec", help="A space separated list of species labels that make up the polyploid clade. Example: 'x,y,z y,z'", default=False);
	parser.add_argument("-h2", dest="h2_spec", help="A space separated list of species labels that make up the copy clade. If spec tree type (-t) is m, this option can be ignored. Example: 'c'", default=False);
	parser.add_argument("-c", dest="group_cap", help="The maxmimum number of groups to consider for any gene tree. Default: 8. Max value: 15.", type=int, default=8);
	parser.add_argument("-o", dest="output_file", help="Output file name.")
	parser.add_argument("-v", dest="verbosity", help="An option to control the amount of output printed to the screen. 0: print only a progress bar. 1: print some output. Default: 1", type=int, default=1);
	parser.add_argument("--checknum", dest="check_num", help="Use this flag in conjunction with all other options to check the number of nodes, groups, and combinations for each gene tree and MUL-tree. In general, gene trees with more than 15 groups to map take a very long time.", action="store_true");
	parser.add_argument("--labeltree", dest="label_opt", help="If this flag is set, the program will read your species tree and simply print it out with the internal nodes labeled.", action="store_true");
	parser.add_argument("--multree", dest="mul_opt", help="Use this along with -s and possibly -h1 and -h2 to simply build MUL-trees from those options.", action="store_true");
	parser.add_argument("--tests", dest="test_opt", help="Use 'grampa.py --tests' the first time you run grampa to run through all the options with pre-set input files.", action="store_true");
	parser.add_argument("--simpson", dest="s_opt", help=argparse.SUPPRESS, action="store_true");

	args = parser.parse_args();

	if args.s_opt:
		RC.simpson();
		sys.exit();

	if args.test_opt:
		t_path = os.path.join(os.path.dirname(__file__), "lib", "tests.py");
		os.system("python " + t_path + " " + sys.version[:3]);
		sys.exit();

	if errorflag == 0:
		if args.label_opt and args.spec_tree == None:
			RC.errorOut(1, "When --labeltree is set, -s must also be set");
			optParse(1);

		elif args.label_opt:
			if args.verbosity != -1:
				print("\n*** Message: --labeltree is set to True! Just labeling your species (-s) tree. All other options will be ignored!");
		elif args.mul_opt and args.spec_tree == None:
			RC.errorOut(2, "When --multree is set, -s must also be set");
			optParse(1);

		elif args.mul_opt:
			if args.verbosity != -1:
				print("\n*** Mesage: --multree is set to True! Just printing out all the MUL-trees you requested. All options but -s, -h1, and -h2 will be ignored!");

		else:
			if args.spec_tree == None or args.spec_tree_type == None or args.gene_input == None or args.h1_spec == None:
				RC.errorOut(3, "-s and -g must be specified.");
				optParse(1);

			if args.output_file == None:
				args.output_file = "MTR_out_" + RC.getLogTime() + ".txt";

			if args.spec_tree_type.lower() not in ['m', 's']:
				RC.errorOut(4, "-m must take values of either m or s");
				optParse(1);

			if args.spec_tree_type.lower() == "m" and (args.h1_spec != False or args.h2_spec != False):
				print("*** Message: With a MUL-tree as the input species tree (-t m) input for -h1 and -h2 are not required.");
				print("*** Your input for -h1 and -h2 will be ignored.")

			if args.group_cap > 15:
				RC.errorOut(5, "For computational reasons, -p should not be set higher than 15.");
				optParse(1);
			elif args.group_cap >=10:
				print("*** Warning! With -p set to 10 or higher, some gene trees may take a very long time to reconcile!");

			if args.verbosity not in [0,1,-1,-2]:
				RC.errorOut(6, "-v must take values of either 0, or 1");
				optParse(1);

		return args.spec_tree, args.spec_tree_type.lower(), args.gene_input, args.h1_spec, args.h2_spec, args.group_cap, args.output_file, args.verbosity, args.check_num, args.label_opt, args.mul_opt;

	elif errorflag == 1:
		parser.print_help();
		print()
		sys.exit();

############################################
def hInParse(h_list):
# Parses the input h nodes.

	if h_list:
		if " " in h_list:
			h_list = h_list.split(" ");
			h_list = list(map(set, [tmp_h.split(",") for tmp_h in h_list]));
		else:
			h_list = list(map(set, [h_list.split(",")]));

	return h_list;

############################################
def getHNodes(spec_list, tree_dict):
# Takes input clades or MUL-tree clades and gets the ancestral node or errors out appropriately.

	h_node, h_mono = RT.LCA(spec_list, tree_dict);
	if not h_mono:
		RC.errorOut(15, "All hybrid clades specified (either in your input MUL-tree or with h1 and h2) must be monophyletic!")
		optParse(1);
	return h_node;

############################################
#Main Block
############################################

starttime = time.time();

spec_file, spec_type, gene_file, hybrid_clades, copy_clades, cap, outfilename, v, check_nums, lab_opt, mul_opt = optParse(0);
# Getting the input parameters.

hybrid_clades = hInParse(hybrid_clades);
copy_clades = hInParse(copy_clades);

if v != -1:
	print("# Reading species tree...");
try:
	spec_tree = open(spec_file, "r").read().replace("\n", "").replace("\r","");
	spec_tree = RT.remBranchLength(spec_tree);
	spec_check = spec_tree.replace("(","").replace(")","").replace(";","").split(",");

	if spec_type == "m":
		hybrid_spec = set([tip for tip in spec_check if spec_check.count(tip) != 1]);
		for h in hybrid_spec:
			spec_tree = spec_tree.replace(h, h+"*", 1);
		hybrid_clades = [hybrid_spec];
	# If the user entered a MUL-tree, some internal re-labeling must be done to those labels that appear twice.

	sinfo, st = RT.treeParse(spec_tree);
	# Parsing of the species tree.

except:
	RC.errorOut(7, "Error reading species tree file!");
	optParse(1);
# Reading the species tree file.

if lab_opt:
	if v != -1:
		print("# The input species tree with internal nodes labeled:");
		print(st + "\n");
	sys.exit();
# The output if --labeltree is set.

if RC.hErrorCheck(spec_check, sinfo, spec_type, hybrid_clades, copy_clades) == 1:
	optParse(1);
# Error handling to make sure the species the user entered are all in their species tree.
## End species tree block.

if not mul_opt:
	if v != -1:
		print("# Reading gene trees from file...");
	try:
		gene_trees = open(gene_file, "r").readlines();
	except:
		RC.errorOut(13, "Error reading gene trees file!");
		optParse(1);
# Reading the gene trees file.
## End gene tree block.

if not check_nums and not mul_opt:
	if outfilename[-4:] == ".txt":
		detoutfilename = outfilename.replace(".txt", "_det.txt");
	else:
		detoutfilename = outfilename + "_det.txt";

	detoutfile = open(detoutfilename, "w");
	detoutfile.write("");
	detoutfile.close();
# If --checknum is not set, we have to prepare both the main and detailed output files.

outfile = open(outfilename, "w");
outfile.write("");
outfile.close();
# Output file prep.
## End output file block.

### Begin input info block!
pad = 65

main_v = 1;
if v == -1:
	main_v = 0;

RC.printWrite(outfilename, main_v, "# =========================================================================");
RC.printWrite(outfilename, main_v, "#\t\t\tMUL-tree reconciliation");
RC.printWrite(outfilename, main_v, "#\t\t\t" + RC.getDateTime());
RC.printWrite(outfilename, main_v, "# The input species tree with internal nodes labeled:", st, pad);
RC.printWrite(outfilename, main_v, "# Main results will be written to file:", outfilename, pad);
if not check_nums and not mul_opt:
	RC.printWrite(outfilename, main_v, "# Detailed results will be written to file:", detoutfilename, pad);

## Identifying the hybrid and copy node in the species tree.
hybrid_nodes = [];
copy_nodes = [];
# print hybrid_clades;
# print copy_clades;
if spec_type == 's':
	RC.printWrite(outfilename, main_v, "# Input species tree is:", "Standard", pad);

	if not hybrid_clades:
		hybrid_nodes = list(sinfo.keys());
		RC.printWrite(outfilename, main_v, "# No H1 node defined", "Searching all possible H1 nodes.", pad);
	else:
		for hybrid_clade in hybrid_clades:
			hybrid_clade = list(hybrid_clade);
			if hybrid_clade[0].isdigit():
				hybrid_node = "<" + hybrid_clade[0] + ">";
			else:
				hybrid_node = getHNodes(hybrid_clade, sinfo);

			if hybrid_node not in hybrid_nodes:
				hybrid_nodes.append(hybrid_node);
		RC.printWrite(outfilename, main_v, "# H1 node(s) identified as:", ",".join(hybrid_nodes), pad);

	if not copy_clades:
		copy_nodes = list(sinfo.keys());
		RC.printWrite(outfilename, main_v, "# No H2 node defined", "Searching all possible H2 nodes.", pad);
	else:
		for copy_clade in copy_clades:
			copy_clade = list(copy_clade);
			if copy_clade[0].isdigit():
				copy_node = "<" + copy_clade[0] + ">";
			else:
				copy_node = getHNodes(copy_clade, sinfo);

			if copy_node not in copy_nodes:
				copy_nodes.append(copy_node);
		RC.printWrite(outfilename, main_v, "# H2 node(s) identified as:", ",".join(copy_nodes), pad);

elif spec_type == 'm':
	hybrid_node = getHNodes(list(hybrid_clades[0]), sinfo);
	copy_node = getHNodes([spec + "*" for spec in list(hybrid_clades[0])], sinfo);

	RC.printWrite(outfilename, main_v, "# Input species tree is:", "MUL-tree", pad);
	RC.printWrite(outfilename, main_v, "# H1 node identified as:", hybrid_node, pad);
	RC.printWrite(outfilename, main_v, "# H2 node identified as:", copy_node, pad);

	hybrid_nodes.append(hybrid_node);
	copy_nodes.append("holder");
	# If the input species tree type is set to 'm', then no copy nodes need to be specified. A holder variable is
	# used as a flag so no MUL-tree is built later.

## End hybrid and copy node identification block.

if check_nums:
	RC.printWrite(outfilename, main_v, "# --checknums set. NOT doing reconciliations, just running some numbers for you...");
elif mul_opt:
	RC.printWrite(outfilename, main_v, "# --multree set. NOT doing reconciliations, just building your MUL-trees...");
RC.printWrite(outfilename, main_v, "# ---------");
### End input info block!

if v == 0:
	if not check_nums:
		print("Beginning reconciliations...\n")

	numiters = len(hybrid_nodes) * len(copy_nodes) * len(gene_trees);
	numbars = 0;
	donepercent = [];
	itercount = 0;
	# Stuff for the loading bar...

mul_dict = {};
mul_num = 1;
hybrid_num = 0;
# Keeping track of scores for each MUL-tree

gt_groups = {};

for hybrid_node in hybrid_nodes:
	if sinfo[hybrid_node][2] == 'root':
		continue;

	if v == -2:
		print("h1", hybrid_node);

	hybrid_num += 1;
	#RC.printWrite(outfilename, 0, "H1-" + str(hybrid_num) + " Node:\t" + hybrid_node);
	if not check_nums and not mul_opt:
		RC.printWrite(detoutfilename, 0, "H1-" + str(hybrid_num) + " Node:\t" + hybrid_node);

	hybrid_clade = set(RT.getClade(hybrid_node, sinfo));

	group_flag = 1;

	for copy_node in copy_nodes:
		if v == -2:
			print("h2", copy_node);

		if copy_node != "holder":
		# Build the MUL-tree if the user has entered a non-MUL species tree.
			if v == 1:
				if len(copy_nodes) == 1:
					print("# Building MUL-tree...");
				else:
					print("# Building MUL-tree for copy node: " + copy_node);

			mt_unlabel = RT.buildMultree(hybrid_node, copy_node, st, sinfo);
			# Building the MUL-tree by passing the species tree to the buildMultree function
			# Input is one node at which to copy the subtree (hybrid_node)
			# and one node at which to place the copy (copy_node)

			if mt_unlabel == "NULL":
				if v == 1:
					print("# Copy node within hybrid subtree or at root. Skipping node.");
					print("# ---------------------------");
				continue;
			# Nodes within the hybrid subtree cannot be copy nodes.

			minfo, mt = RT.treeParse(mt_unlabel);
			# Labeling the MUL-tree as usual with RT.
			if mul_opt:
				outline = hybrid_node + "\t" + copy_node + "\t" + RT.mulPrint(mt, hybrid_clade);
				RC.printWrite(outfilename, v, outline);
				continue;

			if v == -2:
				print("# The MUL-tree w/o internal nodes labeled:", mt_unlabel);
				print("# The MUL-tree with internal nodes labeled:", mt);
				print("# ---------------------------");
				sys.exit();

			mul_dict[mul_num] = [mt, hybrid_clade, copy_node, 0];
			# mul_dict stores, for each mul_tree, the tree, the copy node, and the summed mutation score over all gene trees.
			if check_nums:
				RC.printWrite(outfilename, v, "MT-" + str(mul_num) + "\tTree:" + RT.mulPrint(mt, hybrid_clade) + "\tCopyNode:" + copy_node);
			#else:
			#	RC.printWrite(detoutfilename, v, "MT-" + str(mul_num) + "\tTree:" + mt + "\tCopyNode:" + copy_node);

		else:
		# If the user entered a MUL-tree as their species tree, just assign it here.
			minfo = sinfo;
			mt = st;
			mul_dict[mul_num] = [mt, "", "", 0];

		gene_num = 0;
		num_skipped = 0;

		if check_nums:
			RC.printWrite(outfilename, v, "# Groups\t# Fixed\t# Combinations")

		for gene_tree in gene_trees:
			if v == 0:
				numbars, donepercent = RC.loadingBar(itercount, numiters, donepercent, numbars);
				itercount = itercount + 1;
			# Only the loading bar displays when the program is running if -v is set to 0.

			gene_num = gene_num + 1;

			if gene_tree.strip() == '':
				if not check_nums:
					RC.printWrite(detoutfilename, v, "GT-" + str(gene_num) + "\tEmpty line -- skipping.");
				num_skipped += 1;
				continue;

			try:
				gene_tree = RT.remBranchLength(gene_tree);
				ginfo, gt = RT.treeParse(gene_tree);
			except:
				if not check_nums:
					RC.printWrite(detoutfilename, v, "GT-" + str(gene_num) + "\tError reading this tree! -- skipping.");
				num_skipped += 1;
				continue;

			if len([g for g in ginfo if ginfo[g][2] != 'tip']) != len([g for g in ginfo if ginfo[g][2] == 'tip']) - 1:
				if not check_nums:
					RC.printWrite(detoutfilename, v, "GT-" + str(gene_num) + "\tThis line may not contain a tree, or if so it may be unrooted -- skipping.");
				else:
					print("GT-" + str(gene_num) + "\tThis line may not contain a tree, or if so it may be unrooted -- skipping.");
				num_skipped += 1;
				continue;
			# Parsing the current gene tree.

			if group_flag == 1:
				cur_groups = ALG.collapseGroups(ginfo, hybrid_clade, v);
				gt_groups[gene_num] = cur_groups;
				#print cur_groups;

			outline = "GT-" + str(gene_num) + " to MT-" + str(mul_num) + "\t";
			# Parsing the gene tree.
			if v == -2:
				print('gene tree:', gene_tree);
				print('gt:', gt);
				print("ginfo:", ginfo);
				#sys.exit();
			if not check_nums:
				dup_score, loss_score, maps = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, gt_groups[gene_num], cap, v, check_nums);
			else:
				num_groups, num_fixed = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, gt_groups[gene_num], cap, v, check_nums);
				outline += str(num_groups) + "\t" + str(num_fixed) + "\t" + str(2**num_groups);
				if num_groups > cap:
					outline += " -- Over Cap";
				RC.printWrite(outfilename, v, outline);
				continue;
			# The call of the reconciliation algorithm! On the current gene tree with the current MUL-tree.

			if dup_score == "OVER":
				outline += "Number of groups (" + str(loss_score) + ") over group cap (-p set to " + str(cap) + ") -- Skipping."
				num_skipped += 1;
			else:
				mut_score = dup_score + loss_score;
				outline = outline + str(dup_score) + "\t" + str(loss_score) + "\t" + str(mut_score);
				mul_dict[mul_num][3] += mut_score;
			RC.printWrite(detoutfilename, v, outline);
			# Output stats for the current reconciliation.

		if not check_nums and not mul_opt:
			RC.printWrite(detoutfilename, v, "Number of trees skipped: " + str(num_skipped));
			RC.printWrite(detoutfilename, v, "Total parsimony score for MT-" + str(mul_num) + ": " + str(mul_dict[mul_num][3]));
			RC.printWrite(detoutfilename, v, "# ---------------------------");

			RC.printWrite(outfilename, 0, "MT-" + str(mul_num) + "\t" + hybrid_node + "\t" + copy_node + "\t" + RT.mulPrint(mt, hybrid_clade) + "\t" + str(mul_dict[mul_num][3]));

		# Output total score for the current MUL-tree

		mul_num = mul_num + 1;
		group_flag = 0;

	if not check_nums and not mul_opt:
		min_score = 999999;
		min_tree = "";
		min_num = 0;

		for mtree in mul_dict:
			if mul_dict[mtree][3] < min_score:
				min_score = mul_dict[mtree][3];
				min_tree = mul_dict[mtree][0];
				min_clade = mul_dict[mtree][1];
				min_num = mtree;
		# Look through all the MUL-trees and find the one with the minimum score.

		if v == 0:
			pstring = "100.0% complete.";
			sys.stderr.write('\b' * len(pstring) + pstring);

		RC.printWrite(detoutfilename, 0, "The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + RT.mulPrint(min_tree, min_clade));
		RC.printWrite(outfilename, 0, "# ---------")		
		RC.printWrite(detoutfilename, 0, "Score = " + str(min_score));
		
if not check_nums and not mul_opt:
	RC.printWrite(outfilename, 0, "# The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + RT.mulPrint(min_tree, min_clade));
	RC.printWrite(outfilename, 0, "# Score = " + str(min_score));
	RC.printWrite(detoutfilename, main_v, "\nDone!\n");
endtime = time.time();
totaltime = endtime - starttime;
RC.printWrite(outfilename, main_v, "# Total execution time: " + str(round(totaltime,3)) + " seconds.");
RC.printWrite(outfilename, main_v, "# =========================================================================");





