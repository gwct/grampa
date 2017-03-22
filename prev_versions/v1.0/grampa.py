#!/usr/bin/python
#############################################################################
# Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.
# This is the main interface and handles user options and tree searching.
#
# Gregg Thomas
# Fall 2015
# Combo algorithm implemented Spring 2016
# Filter update implemented Feb. 2017
#############################################################################

import sys, os, re, time, lib.recontree as RT, lib.reconcore as RC, lib.mul_recon as ALG, lib.lca_check as LCHECK

############################################
# Function Definitions
############################################
def optParse(errorflag):
# This function handles the command line options.

	try:
		import argparse;
	except:
		RC.errorOut(0, "Your installation of Python is missing the argparse module. Please try a different version of Python (2.7 or later), or install the module.")
		sys.exit();
	# First check if the argparse module is installed. If not, the input options cannot be parsed.

	parser = argparse.ArgumentParser(description="GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.");

	parser.add_argument("-s", dest="spec_tree", help="A file or string containing a ROOTED, bifurcating, newick formatted species tree in newick format on which to search for polyploid events.");
	parser.add_argument("-t", dest="spec_tree_type", help="[m or s] -- m: input species tree is a MUL-tree. s: input species tree is a singly-labeled tree. Default: s", default="s");
	parser.add_argument("-g", dest="gene_input", help="A file containing one or more ROOTED, bifurcating, newick formatted gene trees to reconcile. The labels in the gene tree MUST end with '_[species name]' and contain no other underscores.");
	parser.add_argument("-h1", dest="h1_spec", help="A space separated list of species labels that make up the polyploid clade. Example: 'x,y,z y,z' Or a space separated list of nodes that make up the polyploid clade. Example '2 4'", default=False);
	parser.add_argument("-h2", dest="h2_spec", help="A space separated list of species labels or internal node labels that make up the clade that you wish to place the second polyploid clade sister to. If spec tree type (-t) is m, this option can be ignored. Example: 'c'", default=False);
	parser.add_argument("-c", dest="group_cap", help="The maxmimum number of groups to consider for any gene tree. Default: 8. Max value: 15.", type=int, default=8);
	parser.add_argument("-o", dest="output_file", help="Output file name.")
	parser.add_argument("-v", dest="verbosity", help="An option to control the amount of output printed to the screen. 0: print only a progress bar and some info. 1: print some detailed output for each reconciliation (this detailed output is also available by default in the _det output file). Default: 1", type=int, default=1);
	parser.add_argument("--standard", dest="lca_standard", help="Select this option to only do standard LCA mapping on your singly-labeled species tree.", action="store_true");
	parser.add_argument("--checknum", dest="check_num", help="Use this flag in conjunction with all other options to check the number of nodes, groups, and combinations for each gene tree and MUL-tree. In general, gene trees with more than 15 groups to map take a very long time.", action="store_true");
	parser.add_argument("--labeltree", dest="label_opt", help="If this flag is set, the program will read your species tree and simply print it out with the internal nodes labeled.", action="store_true");
	parser.add_argument("--multree", dest="mul_opt", help="Use this along with -s and possibly -h1 and -h2 to simply build MUL-trees from those options.", action="store_true");
	parser.add_argument("--tests", dest="test_opt", help="Use 'grampa.py --tests' the first time you run grampa to run through all the options with pre-set input files.", action="store_true");
	parser.add_argument("--simpson", dest="s_opt", help=argparse.SUPPRESS, action="store_true");

	args = parser.parse_args();

	if errorflag == 0:
		if args.s_opt:
			RC.simpson();
			sys.exit();
		# ???

		if args.test_opt:
			t_path = os.path.join(os.path.dirname(__file__), "lib", "tests.py");
			pyver = sys.version[:3];
			os.system("python" + pyver + " " + t_path + " " + pyver);
			sys.exit();
		# Call of the tests script if --tests is set.

		if args.output_file == None:
			args.output_file = "GRAMPA_out_" + RC.getLogTime() + ".txt";
		# Sets a default output file name if none is specified.

		if args.label_opt and args.spec_tree == None:
			RC.errorOut(1, "When --labeltree is set, -s must also be set");
			optParse(1);
		# Checks to see if -s is specified for --labeltree.

		elif args.label_opt:
			if args.verbosity != -1:
				print("\n*** Message: --labeltree is set to True! Just labeling your species (-s) tree. All other options will be ignored!\n");

		elif args.lca_standard:
			print("\n*** Message: --standard is set to True! Just doing standard LCA mapping on your singly-labeled species tree. All options but -s, -g, and -o will be ignored!\n");
			if args.spec_tree_type != 's':
				RC.errorOut(14, "-t must be set to s with the --standard flag!");
				optParse(1);

		elif args.mul_opt and args.spec_tree == None:
			RC.errorOut(2, "When --multree is set, -s must also be set");
			optParse(1);

		elif args.mul_opt:
			if args.verbosity != -1:
				print("\n*** Message: --multree is set to True! Just printing out all the MUL-trees you requested. All options but -s, -h1, and -h2 will be ignored!\n");

		else:
			if args.spec_tree == None or args.spec_tree_type == None or args.gene_input == None or args.h1_spec == None:
				RC.errorOut(3, "-s and -g must be specified.");
				optParse(1);

			if args.spec_tree_type.lower() not in ['m', 's']:
				RC.errorOut(4, "-t must take values of either m or s");
				optParse(1);

			if args.spec_tree_type.lower() == "m" and (args.h1_spec != False or args.h2_spec != False):
				print("*** Message: With a MUL-tree as the input species tree (-t m) input for -h1 and -h2 are not required.");
				print("*** Your input for -h1 and -h2 will be ignored.\n")

			if args.group_cap > 18:
				RC.errorOut(5, "For computational reasons, -p should not be set higher than 15.");
				optParse(1);
			elif args.group_cap >=10:
				print("*** Warning! With -p set to 10 or higher, some gene trees may take a very long time to reconcile!\n");

			if args.verbosity not in [0,1,-1,-2]:
				RC.errorOut(6, "-v must take values of either 0, or 1");
				optParse(1);

		return args.spec_tree, args.spec_tree_type.lower(), args.gene_input, args.h1_spec, args.h2_spec, args.group_cap, args.output_file, args.verbosity, args.lca_standard, args.check_num, args.label_opt, args.mul_opt;

	elif errorflag == 1:
		parser.print_help();
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
# Main Block
############################################
starttime = time.time();

spec_file, spec_type, gene_file, hybrid_clades, copy_clades, cap, outfilename, v, lca_opt, check_nums, lab_opt, mul_opt = optParse(0);
# Getting the input parameters.

hybrid_clades = hInParse(hybrid_clades);
copy_clades = hInParse(copy_clades);
# Parses the input h1 and h2 clades.

### Begin input reading block
if v != -1:
	print("# Reading species tree...");
try:
	if os.path.isfile(spec_file):
		spec_tree = open(spec_file, "r").read().replace("\n", "").replace("\r","");
	else:
		spec_tree = spec_file;
	# If the input string is a filename, read the file. Otherwise, just try it as a newick string.

	spec_tree = RT.remBranchLength(spec_tree);
	spec_check = spec_tree.replace("(","").replace(")","").replace(";","").split(",");

	if spec_type == "m":
		hybrid_spec = set([tip for tip in spec_check if spec_check.count(tip) != 1]);
		for h in hybrid_spec:
			spec_tree = spec_tree.replace(h, h+"*", 1);
		hybrid_clades = [hybrid_spec];
		print hybrid_clades;
	# If the user entered a MUL-tree, some internal re-labeling must be done to those labels that appear twice.

	sinfo, st = RT.treeParse(spec_tree);
	# Parsing of the species tree.

except:
	RC.errorOut(7, "Error reading species tree!");
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
outfile = open(outfilename, "w");
outfile.write("");
outfile.close();
outdir = os.path.dirname(os.path.abspath(outfilename));

if not mul_opt:
	if v != -1:
		print("# Reading gene trees from file...");
	try:
		gene_trees = open(gene_file, "r").readlines();
	except:
		RC.errorOut(13, "Error reading gene trees file!");
		optParse(1);

	gene_file_prefix, gene_file_ext = os.path.splitext(gene_file);
	gene_file_filtered = os.path.join(outdir, gene_file_prefix + "_filtered" + gene_file_ext);
	
# Reading the gene trees file.
## End gene tree block.
### End input reading block.
### Begin output file prep.
if not check_nums and not mul_opt:
	out_file_prefix, out_file_ext = os.path.splitext(outfilename);
	detoutfilename = out_file_prefix + "_det" + out_file_ext;
	detoutfile = open(detoutfilename, "w");
	detoutfile.write("");
	detoutfile.close();
# If --checknum is not set, we have to prepare both the main and detailed output files.

if not mul_opt:
	out_file_prefix, out_file_ext = os.path.splitext(outfilename);
	checkfilename = out_file_prefix + "_checknums" + out_file_ext;
	checkfile = open(checkfilename, "w");
	checkfile.write("Trees\t# Groups\t# Fixed\t# Combinations\n");
	checkfile.close();
# Output file prep.
### End output file block.

### Begin input info block!
pad = 65
main_v = 1;

RC.printWrite(outfilename, main_v, "# =========================================================================");
RC.printWrite(outfilename, main_v, "#\t\t\tMUL-tree reconciliation");
RC.printWrite(outfilename, main_v, "#\t\t\t" + RC.getDateTime());
RC.printWrite(outfilename, main_v, "# The input species tree with internal nodes labeled:", st, pad);
RC.printWrite(outfilename, main_v, "# Main results will be written to file:", outfilename, pad);
if not mul_opt:
	RC.printWrite(outfilename, main_v, "# The number of groups for each tree will be calculated:", checkfilename, pad);
	RC.printWrite(outfilename, main_v, "# Filtered trees will be saved (if necessary):", gene_file_filtered, pad);
if not check_nums and not mul_opt:
	RC.printWrite(outfilename, main_v, "# Detailed results will be written to file:", detoutfilename, pad);

## Identifying the hybrid and copy node in the species tree.
hybrid_nodes = [];
copy_nodes = [];
if spec_type == 's':
	RC.printWrite(outfilename, main_v, "# Input species tree is:", "Standard", pad);

	if lca_opt:
		RC.printWrite(outfilename, main_v, "# --standard set. Only doing standard LCA reconcilition to your singly-labeled gene tree!");

	else:
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
	hybrid_clade = list(hybrid_clades[0]);
	sis_copy_node = getHNodes([spec + "*" for spec in list(hybrid_clades[0])], sinfo);
	copy_node = RT.getSis(sis_copy_node, sinfo);

	RC.printWrite(outfilename, main_v, "# Input species tree is:", "MUL-tree", pad);
	RC.printWrite(outfilename, main_v, "# H1 node identified as:", hybrid_node, pad);
	RC.printWrite(outfilename, main_v, "# H2 node identified as:", copy_node, pad);
## End hybrid and copy node identification block.

if check_nums:
	RC.printWrite(outfilename, main_v, "# --checknums set. NOT doing reconciliations, just running some numbers for you...");
elif mul_opt:
	RC.printWrite(outfilename, main_v, "# --multree set. NOT doing reconciliations, just building your MUL-trees...");
RC.printWrite(outfilename, main_v, "# ---------");
### End input info block!
### Building the MUL-trees, unless only standard recon is being done.
mul_dict = {};
if not lca_opt:
	print("# Building " + str(len(hybrid_nodes) * len(copy_nodes)) + " MUL-trees...");

	mul_num = 1;
	gt_groups = {};

	if spec_type == 's':
		for hybrid_node in hybrid_nodes:
			#if sinfo[hybrid_node][2] == 'root':
			#	continue;	
			if v == -2:
				print("h1", hybrid_node);

			hybrid_clade = set(RT.getClade(hybrid_node, sinfo));

			for copy_node in copy_nodes:
				mt_unlabel = RT.buildMultree(hybrid_node, copy_node, st, sinfo);		
				# Building the MUL-tree by passing the species tree to the buildMultree function
				# Input is one node at which to copy the subtree (hybrid_node)
				# and one node at which to place the copy (copy_node)

				if mt_unlabel == "NULL":
					#print("\n*** Warning: H2 node (" + copy_node + ") within hybrid subtree rooted at H1 node (" + hybrid_node + ") or at root of tree. Not building MUL-tree for this combination.\n");
					#print("# ---------------------------");
					continue;
				# Nodes within the hybrid subtree cannot be copy nodes.

				minfo, mt = RT.treeParse(mt_unlabel);
				# Labeling the MUL-tree as usual with RT.
				if mul_opt:
					outline = hybrid_node + "\t" + copy_node + "\t" + RT.mulPrint(mt, hybrid_clade);
					RC.printWrite(outfilename, v, outline);
					continue;

				mul_dict[mul_num] = [mt, minfo, hybrid_clade, hybrid_node, copy_node, 0, []];
				mul_num += 1;
				# mul_dict stores, for each mul_tree, the tree, the copy node, and the summed mutation score over all gene trees.				
	else:
	# If the user entered a MUL-tree as their species tree, just assign it here.
		mul_dict[mul_num] = [st, sinfo, hybrid_clade, hybrid_node, copy_node, 0, []];
else:
	mul_dict["st"] = "standard";
# This block builds the MUL-trees and prepares the main mul_dict:
# mul_dict -> mul_num : [MUL-tree string, MUL-tree dict, hybrid clade species, hybrid node, copy node, an initial score of 0, an empty list to add gene tree groups to later]

if mul_opt:
	endtime = time.time();
	totaltime = endtime - starttime;
	RC.printWrite(outfilename, main_v, "# Total execution time: " + str(round(totaltime,3)) + " seconds.");
	RC.printWrite(outfilename, main_v, "# =========================================================================");
	sys.exit();
## If --mulopt is set, this just prints out the MUL-trees and exits.
### End MUL-tree block.
### Begin gene tree parsing and filtering block.
print("# Checking gene trees...");
checkfile = open(checkfilename, "a");

gene_num = -1;
num_skipped = 0;
gene_trees_filtered = [];
# gene_trees_filtered will be a list of lists. One list for each gene tree. If the gene tree passes all
# filters, it will be [gene tree string, gene tree dict]. Otherwise it will be [Filter message].

for gene_tree in gene_trees:
	gene_num += 1;

	if gene_tree.strip() == '':
		gene_trees_filtered.append(["# Empty line -- Filtering."]);
		checkfile.write("GT-" + str(gene_num+1) + "\tEmpty line -- Filtering.\n");
		num_skipped += 1;
		continue;
	# Handles empty lines in the gene tree file.

	try:
		gene_tree = RT.remBranchLength(gene_tree);
		ginfo, gt = RT.treeParse(gene_tree);
		gene_trees_filtered.append([gt,ginfo]);
	except:
		gene_trees_filtered.append(["# Error reading this tree! -- Filtering."]);
		checkfile.write("GT-" + str(gene_num+1) + "\tError reading this tree! -- Filtering.\n");
		num_skipped += 1;
		continue;
	# Tries the gene tree parsing code and if anything goes wrong, catches exception and filters the tree.

	if len([g for g in ginfo if ginfo[g][2] != 'tip']) != len([g for g in ginfo if ginfo[g][2] == 'tip']) - 1:
		gene_trees_filtered[gene_num] = ["# This line may not contain a tree, or if so it may be unrooted -- Filtering."]
		checkfile.write("GT-" + str(gene_num+1) + "\tThis line may not contain a tree, or if so it may be unrooted -- Filtering.\n");
		num_skipped += 1;
		continue;
	# Another check for gene tree parsing and formatting errors.
## First read through of gene trees to check for parsing errors.

if not lca_opt:
# Do not need to check groups if only doing standard LCA.
	mul_dict, gene_trees_filtered, num_skipped = ALG.collapseGroups(mul_dict, sinfo, gene_trees_filtered, checkfile, num_skipped, cap, v);
	# The call of the important collapseGroups function that groups polyploid clades in the gene trees to speed up the reconciliations to MUL-trees.
	checkfile.close();

	if num_skipped != 0:
		print("# Filtered " + str(num_skipped) + " trees.")
		print("# Writing filtered gene trees to file...");
		filtered_file = open(gene_file_filtered, "w");
		for gt in gene_trees_filtered:
			filtered_file.write(gt[0] + "\n");
		filtered_file.close();
	else:
		print("# No trees filtered! Using your original set.")
	# Write the filtered trees to a file, or not if no filtering was done.

	if check_nums:
		endtime = time.time();
		totaltime = endtime - starttime;
		RC.printWrite(outfilename, main_v, "# Total execution time: " + str(round(totaltime,3)) + " seconds.");
		RC.printWrite(outfilename, main_v, "# =========================================================================");
		sys.exit();
	# If --checknums is set, exit the program here.
### End gene tree block.
### Begin standard recon block
if spec_type == 's':
# Only do standard recon if the input tree is a singly-labeled tree.
	print("# Doing standard reconciliation on your singly-labeled tree...\n");
	st_score = 0;
	gene_num = -1;

	tot_node_counts = {};
	for node in sinfo:
		tot_node_counts[node] = [0,0];

	RC.printWrite(detoutfilename, v, "# ---------------------------");
	RC.printWrite(detoutfilename, v, "ST\t" + st);
	for gene_tree in gene_trees_filtered:
		gene_num += 1;
		if len(gene_tree) == 1:
			continue;
		# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

		outline = "GT-" + str(gene_num+1) + " to ST\t";
		gt, ginfo = gene_tree;
		# Retrieves the gene tree info.

		maps = {};
		for g in ginfo:
			if ginfo[g][2] == 'tip':
				speclabel = g[g.rfind("_")+1:];
				maps[g] = [speclabel];
			else:
				maps[g] = [];
		# Initialize the maps.

		st_maps, st_num_dups, st_num_loss, st_node_counts = ALG.reconLCA(ginfo, sinfo, maps);
		st_mut_score = st_num_dups + st_num_loss;
		st_score += st_mut_score;
		outline = outline + str(st_num_dups) + "\t" + str(st_num_loss) + "\t" + str(st_mut_score);
		RC.printWrite(detoutfilename, v, outline);
		# Call the recon algorithm and aggregate scores.

		# branch_outline = "\t";
		# for node in sorted(st_node_counts.keys()):
		# 	tot_node_counts[node][0] += st_node_counts[node][0];
		# 	tot_node_counts[node][1] += st_node_counts[node][1];

		# 	branch_outline += "\t" + node + ":" + str(st_node_counts[node][0]) + "," + str(st_node_counts[node][1])
		# RC.printWrite(detoutfilename, v, branch_outline);
		# Write the branch gain/loss scores.

	RC.printWrite(detoutfilename, v, "Total parsimony score for ST: " + str(st_score));
	# branch_outline = "Total branch scores for ST:\t" + "\t".join([node + ":" + str(tot_node_counts[node][0]) + "," + str(tot_node_counts[node][1]) for node in sorted(tot_node_counts.keys())]);
	# RC.printWrite(detoutfilename, v, branch_outline);
	RC.printWrite(detoutfilename, v, "# ---------------------------");
	RC.printWrite(outfilename, 0, "ST\t\t\t" + st + "\t" + str(st_score));
 	# Print the total score and total branch scores for the singly-labeled tree.

if lca_opt:
	print("# Done!");
	sys.exit();
# If --standard is set, exit the program here.
### End standard recon block.
### Begin MUL-recon block.
print("# Beginning reconciliations on filtered gene trees...\n")
if v == 0:
	numiters = len(mul_dict) * len(gene_trees_filtered);
	numbars = 0;
	donepercent = [];
	itercount = 0;
# Stuff for the loading bar...

for mul_num in mul_dict:
	mt = mul_dict[mul_num][0];
	minfo = mul_dict[mul_num][1];
	hybrid_clade = mul_dict[mul_num][2];
	hybrid_node = mul_dict[mul_num][3];
	copy_node = mul_dict[mul_num][4];
	group_list = gt_groups = mul_dict[mul_num][6];

	tot_node_counts = {};
	for node in minfo:
		tot_node_counts[node] = [0,0];

	RC.printWrite(detoutfilename, 0, "MT-" + str(mul_num) + ":" + RT.mulPrint(mt, hybrid_clade) + "\tH1 Node:" + hybrid_node + "\tH2 Node:" + copy_node);
	gene_num = -1;

	for gene_tree in gene_trees_filtered:
		if v == 0:
			numbars, donepercent = RC.loadingBar(itercount, numiters, donepercent, numbars);
			itercount = itercount + 1;
		# Only the loading bar displays when the program is running if -v is set to 0.

		gene_num += 1;
		if len(gene_tree) == 1:
			continue;
		# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

		gt_groups = group_list[gene_num][0];
		gt_fixed = group_list[gene_num][1];
		outline = "GT-" + str(gene_num+1) + " to MT-" + str(mul_num) + "\t";
		gt, ginfo = gene_tree;
		# Retrieve gene tree info and collapsed groups for this gene tree-MUL-tree combo

		dup_score, loss_score, maps, mt_node_counts = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, gt_groups, gt_fixed, cap, v, check_nums);
		# The call of the reconciliation algorithm! On the current gene tree with the current MUL-tree.

		mut_score = dup_score + loss_score;
		outline = outline + str(dup_score) + "\t" + str(loss_score) + "\t" + str(mut_score);
		mul_dict[mul_num][5] += mut_score;
		RC.printWrite(detoutfilename, v, outline);
		# Aggregate scores.

		# branch_outline = "\t";
		# for node in sorted(mt_node_counts.keys()):
		# 	tot_node_counts[node][0] += mt_node_counts[node][0];
		# 	tot_node_counts[node][1] += mt_node_counts[node][1];

		# 	branch_outline += "\t" + node + ":" + str(mt_node_counts[node][0]) + "," + str(mt_node_counts[node][1])
		# RC.printWrite(detoutfilename, v, branch_outline);
		# Write the branch gain/loss scores.

	RC.printWrite(detoutfilename, v, "Total parsimony score for MT-" + str(mul_num) + ": " + str(mul_dict[mul_num][5]));
	# branch_outline = "Total branch scores for MT" + str(mul_num) + ":\t" + "\t".join([node + ":" + str(tot_node_counts[node][0]) + "," + str(tot_node_counts[node][1]) for node in sorted(tot_node_counts.keys())]);
	# RC.printWrite(detoutfilename, v, branch_outline);
	RC.printWrite(detoutfilename, v, "# ---------------------------");
	RC.printWrite(outfilename, 0, "MT-" + str(mul_num) + "\t" + hybrid_node + "\t" + copy_node + "\t" + RT.mulPrint(mt, hybrid_clade) + "\t" + str(mul_dict[mul_num][5]));
	# Print the total score and total branch scores for the current MUL-tree.
### End MUL-recon block.
### Begin scoring block.
min_score = 999999;
min_tree = "";
min_num = 0;

for mul_num in mul_dict:
	if mul_dict[mul_num][5] < min_score:
		min_score = mul_dict[mul_num][5];
		min_tree = mul_dict[mul_num][0];
		min_clade = mul_dict[mul_num][2];
		min_num = mul_num;
# Look through all the MUL-trees and find the one with the minimum score.

if spec_type == 's' and st_score < min_score:
	RC.printWrite(outfilename, 0, "# The tree with the minimum parsimony score is the singly-labled tree (ST):\t" + st);
	RC.printWrite(outfilename, 0, "# Score = " + str(st_score));
else:
	RC.printWrite(outfilename, 0, "# The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + RT.mulPrint(min_tree, min_clade));
	RC.printWrite(outfilename, 0, "# Score = " + str(min_score));
if v == 0:
	pstring = "100.0% complete.";
	sys.stderr.write('\b' * len(pstring) + pstring);
RC.printWrite(detoutfilename, main_v, "\n# Done!\n");
# Final output block.

endtime = time.time();
totaltime = endtime - starttime;
RC.printWrite(outfilename, main_v, "# Total execution time: " + str(round(totaltime,3)) + " seconds.");
RC.printWrite(outfilename, main_v, "# =========================================================================");
### End scoring block and end program!












