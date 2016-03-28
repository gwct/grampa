#!/usr/bin/python
#############################################################################
# Gene tree reconciliation using multi-labeled trees (MUL-trees) to account
# for polyploid events. This interface handles user options and tree searching.
#
# Gregg Thomas
# Fall 2015, Combo algorithm implemented Spring 2016
#############################################################################

import sys, argparse, lib.recontree as RT, lib.reconcore as RC, lib.mul_recon as ALG

# h "51914,CANGA,KAZAF,588726,51660,SACCA,NAUDC,SACBA,YEAST,1071379,TETPH,VANPO"
# c "LACTH,KLUWA,SACKL,ASHGO,KLULA"
# python .\interface.py -s .\trees\s.tre -y "x,y,z" -o .\testout.txt -g .\trees\g05.tre -c "b"

############################################
#Function Definitions
############################################
def optParse(errorflag):
# This function handles the command line options.

	parser = argparse.ArgumentParser(description="Tree reconciliation using MUL-trees to find polyploidy events.");

	parser.add_argument("-s", dest="spec_tree", help="A bifurcating species tree in newick format on which to search for polyploid events.");
	parser.add_argument("-t", dest="spec_tree_type", help="[m or s] -- m: input species tree is a MUL-tree. s: input species tree is a non-MUL tree.")
	parser.add_argument("-g", dest="gene_input", help="A file containing one or more newick formatted gene trees to reconcile. The labels in the gene tree must end with '_[species name]' and contain no other underscores.");
	parser.add_argument("-y", dest="hybrid_spec", help="A comma separated list of species labels that make up the polyploid clade. Example: 'x,y,z'");
	parser.add_argument("-c", dest="copy_spec", help="A comma separated list of species labels that make up the copy clade. If spec tree type (-t) is m, this option can be ignored. Example: 'c'", default="");
	parser.add_argument("-o", dest="output_file", help="Output file name.")
	parser.add_argument("-v", dest="verbosity", help="An option to control the amount of output printed to the screen. 0: print only a progress bar. 1: print some output. Default: 1", type=int, default=1);
	parser.add_argument("--checknum", dest="check_num", help="Use this flag in conjunction with all other options to check the number of nodes, groups, and combinations for each gene tree and MUL-tree. In general, gene trees with more than 20 groups to map take a very long time.", action="store_true");

	args = parser.parse_args();

	if errorflag == 0:
		if args.spec_tree == None or args.spec_tree_type == None or args.gene_input == None or args.hybrid_spec == None:
			RC.errorOut(1, "-s, -t, -g, and -y must all be specified.");
			optParse(1);

		if args.output_file == None:
			args.output_file = "MTR_out_" + RC.getLogTime() + ".txt";

		if args.spec_tree_type.lower() not in ['m', 's']:
			RC.errorOut(2, "-m must take values of either m or s");
			optParse(1);

		if args.spec_tree_type.lower() == "m" and args.copy_spec != None:
			print "*** Message: With a MUL-tree as the input species tree (-t m) a copy node (-c) is not required.";
			print "*** Your input for -c will be ignored."

		if args.verbosity not in [0,1,-2]:
			RC.errorOut(3, "-v must take values of either 0, or 1");
			optParse(1);

		return args.spec_tree, args.spec_tree_type.lower(), args.gene_input, args.hybrid_spec.replace(" ","").split(","), args.copy_spec.replace(" ","").split(","), args.output_file, args.verbosity, args.check_num;

	elif errorflag == 1:
		parser.print_help();
		print
		sys.exit();

############################################
#Main Block
############################################

spec_file, spec_type, gene_file, hybrid_list, copy_list, outfilename, v, check_nums = optParse(0);
# Getting the input parameters.
print check_nums;
print "# Reading species tree...";
try:
	spec_tree = open(spec_file, "r").read().replace("\n", "").replace("\r","");
	spec_tree = RT.remBranchLength(spec_tree);

	if spec_type == "m":
		for h in hybrid_list:
			spec_tree = spec_tree.replace(h, h+"*", 1);
	# If the user entered a MUL-tree, some internal re-labeling must be done to those labels that appear twice.

	sinfo, st = RT.treeParseNew(spec_tree,2);
	# Parsing of the species tree.
except:
	RC.errorOut(4, "Error reading species tree file!");
	optParse(1);
# Reading the species tree file.

### Begin error handling block.
spec_check = spec_tree.replace("(","").replace(")","").split(",");
if not all(h in spec_check for h in hybrid_list):
	RC.errorOut(5, "Not all hybrid species (-y) are present in your species tree!");
	optParse(1);

if copy_list != [''] and not all(c in spec_check for c in copy_list):
	RC.errorOut(6, "Not all copy species (-c) are present in your species tree!");
	optParse(1);

if spec_type == 's' and any(spec_check.count(n) > 1 for n in spec_check):
	RC.errorOut(7, "You have entered a tree type (-t) of 's' but there are labels in your tree that appear more than once!");
	optParse(1);

if spec_type == 'm' and any(spec_check.count(h) != 2 for h in hybrid_list):
	RC.errorOut(8, "You have entered a tree type (-t) of 'm', so all your hybrid species (-h) should appear exactly twice.");
	optParse(1);

if spec_type == 'm' and any(spec_check.count(s) != 1 for s in spec_check if s not in hybrid_list):
	RC.errorOut(9, "You have entered a tree type (-t) of 'm', so all non-hybrid species should appear exactly once.");
	optParse(1);

# Error handling to make sure the species the user entered are all in their species tree.
### End error handling block.

print "# Reading gene trees from file...";
try:
	gene_trees = open(gene_file, "r").readlines();
except:
	RC.errorOut(10, "Error reading gene trees file!");
	optParse(1);
# Reading the gene trees file.

outfile = open(outfilename, "w");
outfile.write("");
outfile.close();
# Output file prep.

### Begin input info block!
hybrid_clade = set(hybrid_list);
pad = 65

RC.printWrite(outfilename, 1, "# =========================================================================");
RC.printWrite(outfilename, 1, "#\t\t\tMUL-tree reconciliation");
RC.printWrite(outfilename, 1, "#\t\t\t" + RC.getDateTime());
RC.printWrite(outfilename, 1, "# The input species tree with internal nodes labeled:", st, pad);
RC.printWrite(outfilename, 1, "# The polyploid clade:", ", ".join(hybrid_list), pad);
RC.printWrite(outfilename, 1, "# Result will be written to file:", outfilename, pad);

for s in sinfo:
	s_clade = RT.getClade(s, sinfo);
	if hybrid_clade == set(s_clade):
		hybrid_node = s;
		RC.printWrite(outfilename, 1, "# Hybrid node identified as:", hybrid_node, pad);
	if copy_list != [] and set(copy_list) == set(s_clade) and spec_type != 'm':
		copy_nodes = s;
		RC.printWrite(outfilename, 1, "# Copy node identified as:", copy_nodes, pad);
		copy_nodes = [copy_nodes];

if spec_type == 'm':
	copy_nodes = ["holder"];
# If the input species tree type is set to 'm', then no copy nodes need to be specified. A holder variable is
# used as a flag so no MUL-tree is built later.

if copy_list == ['']:
	copy_nodes = sinfo.keys();
	RC.printWrite(outfilename, 1, "# No copy node defined\n# Searching all possible copy nodes.");
# Identifying the hybrid and copy node in the species tree.
if check_nums:
	RC.printWrite(outfilename, 1, "# --checknums set. NOT doing reconciliations, just running some numbers for you...");
RC.printWrite(outfilename, 1, "# ---------");
### End input info block!

if v == 0:
	print "Beginning reconciliations...\n"

	numiters = len(copy_nodes) * len(gene_trees);
	numbars = 0;
	donepercent = [];
	itercount = 0;
	# Stuff for the loading bar...

mul_dict = {};
mul_num = 1;
# Keeping track of scores for each MUL-tree

for copy_node in copy_nodes:
	if copy_node != "holder":
	# Build the MUL-tree if the user has entered a non-MUL species tree.
		if v == 1:
			if len(copy_nodes) == 1:
				print "# Building MUL-tree...";
			else:
				print "# Building MUL-tree for copy node:", copy_node;

		mt_unlabel = RT.buildMultree(hybrid_node, copy_node, st, sinfo);
		# Building the MUL-tree by passing the species tree to the buildMultree function
		# Input is one node at which to copy the subtree (hybrid_node)
		# and one node at which to place the copy (copy_node)


		if mt_unlabel == "NULL":
			if v == 1:
				print "# Copy node within hybrid subtree or at root. Skipping node.";
				print "# ---------";
			continue;
		# Nodes within the hybrid subtree cannot be copy nodes.

		minfo, mt = RT.treeParseNew(mt_unlabel,2);
		# Labeling the MUL-tree as usual with RT.
		if v == 1:
			print "# The MUL-tree w/o internal nodes labeled:", mt_unlabel;
			print "# The MUL-tree with internal nodes labeled:", mt;
			print "# ---------";

		mul_dict[mul_num] = [mt, copy_node, 0];
		# mul_dict stores, for each mul_tree, the tree, the copy node, and the summed mutation score over all gene trees.
		RC.printWrite(outfilename, v, "MT-" + str(mul_num) + "\tTree:" + mt + "\tCopyNode:" + copy_node);

	else:
	# If the user entered a MUL-tree as their species tree, just assign it here.
		minfo = sinfo;
		mt = st;
		mul_dict[mul_num] = [mt, "", 0];

	gene_num = 1;

	if check_nums:
		RC.printWrite(outfilename, v, "Trees\t# Singles\t# Groups\t# Fixed\t# Combinations")

	for gene_tree in gene_trees:
		if v == 0:
			numbars, donepercent = RC.loadingBar(itercount, numiters, donepercent, numbars);
			itercount = itercount + 1;
		# Only the loading bar displays when the program is running if -v is set to 0.

		gene_tree = RT.remBranchLength(gene_tree);

		ginfo, gt = RT.treeParseNew(gene_tree,2);
		# Parsing the current gene tree.

		# if v == 1:
		#	print "# Gene tree " + str(gene_num) + " -- Running MUL-reconciliation algorithm...";
		#
		outline = "GT-" + str(gene_num) + " to MT-" + str(mul_num) + "\t";
		# Parsing the gene tree.
		if v == -2:
			print 'gt:', gt;
		gene_num = gene_num + 1;

		if not check_nums:
			dup_score, loss_score, maps = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, v, check_nums);
		else:
			num_singles, num_groups, num_fixed, num_total = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, v, check_nums);
			outline += str(num_singles) + "\t" + str(num_groups) + "\t" + str(num_fixed) + "\t" + str(2**num_total);
			RC.printWrite(outfilename, v, outline);
			continue;
		mut_score = dup_score + loss_score;
		# The call of the reconciliation algorithm! On the current gene tree with the current MUL-tree.

		outline = outline + str(dup_score) + "\t" + str(loss_score) + "\t" + str(mut_score);
		mul_dict[mul_num][2] += mut_score;
		RC.printWrite(outfilename, v, outline);
		if v == -2:
			print "--------_------";
		# Output stats for the current reconciliation.
		if gene_num == 110 and v == -2:
			sys.exit();

	if not check_nums:
		RC.printWrite(outfilename, v, "Total parsimony score for MT-" + str(mul_num) + ": " + str(mul_dict[mul_num][2]));
		RC.printWrite(outfilename, v, "# ---------");
	# Output total score for the current MUL-tree

	mul_num = mul_num + 1;

if not check_nums:
	min_score = 999999;
	min_tree = "";
	min_num = 0;

	for mtree in mul_dict:
		if mul_dict[mtree][2] < min_score:
			min_score = mul_dict[mtree][2];
			min_tree = mul_dict[mtree][0];
			min_num = mtree;
	# Look through all the MUL-trees and find the one with the minimum score.

	if v == 0:
		pstring = "100.0% complete.";
		sys.stderr.write('\b' * len(pstring) + pstring);

	RC.printWrite(outfilename, 1, "\nDone!\n");

	RC.printWrite(outfilename, 1, "The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + min_tree);
	RC.printWrite(outfilename, 1, "Score = " + str(min_score));
RC.printWrite(outfilename, 1, "# =========================================================================");





