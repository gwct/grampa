#!/usr/bin/python
#############################################################################
# Gene tree reconciliation using multi-labeled trees (MUL-trees) to account
# for polyploid events. This interface handles user options and tree searching.
#
# Gregg Thomas
# Fall 2015, Combo algorithm implemented Spring 2016
#############################################################################

import sys, time, argparse, lib.recontree as RT, lib.reconcore as RC, lib.mul_recon as ALG

# Yeast hybrid and copy nodes:
# 	Hybrid:	"51914,CANGA,KAZAF,588726,51660,SACCA,NAUDC,SACBA,YEAST,1071379,TETPH,VANPO"
# 	Copy:	"LACTH,KLUWA,SACKL,ASHGO,KLULA"
#
# Test commands:
#	Manual checknum:		mul_interface.py -s test/in/manual_sims/s.tre -t s -y "x,y,z" -g test/in/manual_sims/tree_list.txt -o test/out/manual_checknum.txt --checknum
#	Manual test:			mul_interface.py -s test/in/manual_sims/s.tre -t s -y "x,y,z" -g test/in/manual_sims/tree_list.txt -o test/out/manual_full_out.txt
#	Yeast 1 tree checknum:	mul_interface.py -s test/in/yeast_spec_abb.tre -t s -y "51914,CANGA,KAZAF,588726,51660,SACCA,NAUDC,SACBA,YEAST,1071379,TETPH,VANPO" -c "LACTH,KLUWA,SACKL,ASHGO,KLULA" -g test/in/yeast206_trees_rfbr.txt -o test/out/yeast206_checknum_1.txt --checknum
#	Yeast 1 tree test:		mul_interface.py -s test/in/yeast_spec_abb.tre -t s -y "51914,CANGA,KAZAF,588726,51660,SACCA,NAUDC,SACBA,YEAST,1071379,TETPH,VANPO" -c "LACTH,KLUWA,SACKL,ASHGO,KLULA" -g test/in/yeast206_trees_rfbr.txt -o test/out/yeast_out/yeast206.8.1.txt
#	Yeast full checknum:	mul_interface.py -s test/in/yeast_spec_abb.tre -t s -y "51914,CANGA,KAZAF,588726,51660,SACCA,NAUDC,SACBA,YEAST,1071379,TETPH,VANPO" -g test/in/yeast206_trees_rfbr.txt -o test/out/yeast206_checknum.txt --checknum
#	Yeast full test:		mul_interface.py -s test/in/yeast_spec_abb.tre -t s -y "51914,CANGA,KAZAF,588726,51660,SACCA,NAUDC,SACBA,YEAST,1071379,TETPH,VANPO" -g test/in/yeast206_trees_rfbr.txt -o test/out/yeast_out/yeast206.8.full.txt

############################################
#Function Definitions
############################################
def optParse(errorflag):
# This function handles the command line options.

	parser = argparse.ArgumentParser(description="Tree reconciliation using MUL-trees to find polyploidy events.");

	parser.add_argument("-s", dest="spec_tree", help="A bifurcating species tree in newick format on which to search for polyploid events.");
	parser.add_argument("-t", dest="spec_tree_type", help="[m or s] -- m: input species tree is a MUL-tree. s: input species tree is a standard tree. Default: s", default="s");
	parser.add_argument("-g", dest="gene_input", help="A file containing one or more newick formatted gene trees to reconcile. The labels in the gene tree must end with '_[species name]' and contain no other underscores.");
	parser.add_argument("-h1", dest="h1_spec", help="A comma separated list of species labels that make up the polyploid clade. Example: 'x,y,z y,z'", default="");
	parser.add_argument("-h2", dest="h2_spec", help="A comma separated list of species labels that make up the copy clade. If spec tree type (-t) is m, this option can be ignored. Example: 'c'", default="");
	parser.add_argument("-p", dest="group_cap", help="The maxmimum number of groups to consider for any gene tree. Default: 8. Max value: 15.", type=int, default=8);
	parser.add_argument("-o", dest="output_file", help="Output file name.")
	parser.add_argument("-v", dest="verbosity", help="An option to control the amount of output printed to the screen. 0: print only a progress bar. 1: print some output. Default: 1", type=int, default=1);
	parser.add_argument("--checknum", dest="check_num", help="Use this flag in conjunction with all other options to check the number of nodes, groups, and combinations for each gene tree and MUL-tree. In general, gene trees with more than 20 groups to map take a very long time.", action="store_true");

	args = parser.parse_args();

	if errorflag == 0:
		if args.spec_tree == None or args.spec_tree_type == None or args.gene_input == None or args.h1_spec == None:
			RC.errorOut(1, "-s and -g must be specified.");
			optParse(1);

		if args.output_file == None:
			args.output_file = "MTR_out_" + RC.getLogTime() + ".txt";

		if args.spec_tree_type.lower() not in ['m', 's']:
			RC.errorOut(2, "-m must take values of either m or s");
			optParse(1);

		if args.spec_tree_type.lower() == "m" and args.h2_spec != None:
			print "*** Message: With a MUL-tree as the input species tree (-t m) a copy node (-c) is not required.";
			print "*** Your input for -c will be ignored."

		if args.group_cap > 15:
			RC.errorOut(3, "For computational reasons, -p should not be set higher than 15.");
			optParse(1);
		elif args.group_cap >=10:
			print "*** Warning! With -p set to 10 or higher, some gene trees may take a very long time to reconcile!";

		if args.verbosity not in [0,1,-2,-3]:
			RC.errorOut(4, "-v must take values of either 0, or 1");
			optParse(1);

		#if " " in args.h1_spec:
		#	args.h1_spec = args.h1_spec.split(" ");
		#ar

		return args.spec_tree, args.spec_tree_type.lower(), args.gene_input, args.h1_spec.replace(" ","").split(","), args.h2_spec.replace(" ","").split(","), args.group_cap, args.output_file, args.verbosity, args.check_num;

	elif errorflag == 1:
		parser.print_help();
		print
		sys.exit();

############################################
#Main Block
############################################

starttime = time.time();

spec_file, spec_type, gene_file, hybrid_list, copy_list, cap, outfilename, v, check_nums = optParse(0);
# Getting the input parameters.

print hybrid_list;


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
	RC.errorOut(5, "Error reading species tree file!");
	optParse(1);
# Reading the species tree file.

### Begin error handling block.
spec_check = spec_tree.replace("(","").replace(")","").replace(";","").split(",");
if hybrid_list != [''] and not all(h in spec_check for h in hybrid_list):
	RC.errorOut(6, "Not all hybrid species (-h1) are present in your species tree!");
	optParse(1);

if copy_list != [''] and not all(c in spec_check for c in copy_list):
	RC.errorOut(7, "Not all copy species (-h2) are present in your species tree!");
	optParse(1);

if spec_type == 's' and any(spec_check.count(n) > 1 for n in spec_check):
	RC.errorOut(8, "You have entered a tree type (-t) of 's' but there are labels in your tree that appear more than once!");
	optParse(1);

if spec_type == 'm' and any(spec_check.count(h) != 2 for h in hybrid_list):
	RC.errorOut(9, "You have entered a tree type (-t) of 'm', so all your hybrid species (-h1) should appear exactly twice.");
	optParse(1);

if spec_type == 'm' and any(spec_check.count(s) != 1 for s in spec_check if s not in hybrid_list):
	RC.errorOut(10, "You have entered a tree type (-t) of 'm', so all non-hybrid species should appear exactly once.");
	optParse(1);

# Error handling to make sure the species the user entered are all in their species tree.
### End error handling block.

print "# Reading gene trees from file...";
try:
	gene_trees = open(gene_file, "r").readlines();
except:
	RC.errorOut(11, "Error reading gene trees file!");
	optParse(1);
# Reading the gene trees file.

if not check_nums:
	if outfilename[-4:] == ".txt":
		detoutfilename = outfilename.replace(".txt", "_det.txt");
	else:
		detoutfilename = outfilename + "_det.txt";

	detoutfile = open(detoutfilename, "w");
	detoutfile.write("");
	detoutfile.close();

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
RC.printWrite(outfilename, 1, "# Main results will be written to file:", outfilename, pad);
if not check_nums:
	RC.printWrite(outfilename, 1, "# Detailed results will be written to file:", detoutfilename, pad);

for s in sinfo:
	s_clade = RT.getClade(s, sinfo);
	if hybrid_list != [''] and hybrid_clade == set(s_clade):
		RC.printWrite(outfilename, 1, "# The polyploid clade:", ", ".join(hybrid_list), pad);
		RC.printWrite(outfilename, 1, "# H1 identified as:", s, pad);
		hybrid_nodes = [s];
	if copy_list != [''] and set(copy_list) == set(s_clade) and spec_type != 'm':
		RC.printWrite(outfilename, 1, "# H2 identified as:", s, pad);
		copy_nodes = [s];

if hybrid_list == ['']:
	hybrid_nodes = sinfo.keys();
	RC.printWrite(outfilename, 1, "# No H1 node defined", "Searching all possible H1 nodes.", pad);
if copy_list == ['']:
	copy_nodes = sinfo.keys();
	RC.printWrite(outfilename, 1, "# No H2 node defined", "Searching all possible H2 nodes.", pad);
# Identifying the hybrid and copy node in the species tree.
if check_nums:
	RC.printWrite(outfilename, 1, "# --checknums set. NOT doing reconciliations, just running some numbers for you...");
RC.printWrite(outfilename, 1, "# ---------");
### End input info block!

if v == 0:
	if not check_nums:
		print "Beginning reconciliations...\n"

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
	itercount = itercount + 1;
	if sinfo[hybrid_node][3] == 'root':
		continue;

	hybrid_num += 1;
	RC.printWrite(outfilename, 0, "H1-" + str(hybrid_num) + " Node:\t" + hybrid_node);
	RC.printWrite(detoutfilename, 0, "H1-" + str(hybrid_num) + " Node:\t" + hybrid_node);

	hybrid_clade = set(RT.getClade(hybrid_node, sinfo));

	group_flag = 1;

	if spec_type == 'm':
		copy_nodes = ["holder"];
	# If the input species tree type is set to 'm', then no copy nodes need to be specified. A holder variable is
	# used as a flag so no MUL-tree is built later.

	for copy_node in copy_nodes:
		itercount = itercount + 1;
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
					print "# ---------------------------";
				continue;
			# Nodes within the hybrid subtree cannot be copy nodes.

			minfo, mt = RT.treeParseNew(mt_unlabel,2);
			# Labeling the MUL-tree as usual with RT.
			if v == -3:
				print "# The MUL-tree w/o internal nodes labeled:", mt_unlabel;
				print "# The MUL-tree with internal nodes labeled:", mt;
				print "# ---------------------------";

			mul_dict[mul_num] = [mt, copy_node, 0];
			# mul_dict stores, for each mul_tree, the tree, the copy node, and the summed mutation score over all gene trees.
			if check_nums:
				RC.printWrite(outfilename, v, "MT-" + str(mul_num) + "\tTree:" + mt + "\tCopyNode:" + copy_node);
			#else:
			#	RC.printWrite(detoutfilename, v, "MT-" + str(mul_num) + "\tTree:" + mt + "\tCopyNode:" + copy_node);

		else:
		# If the user entered a MUL-tree as their species tree, just assign it here.
			minfo = sinfo;
			mt = st;
			mul_dict[mul_num] = [mt, "", 0];

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
				RC.printWrite(detoutfilename, v, "GT-" + str(gene_num) + "\tEmpty line -- skipping.");
				continue;

			gene_tree = RT.remBranchLength(gene_tree);
			ginfo, gt = RT.treeParseNew(gene_tree,2);
			# Parsing the current gene tree.

			if group_flag == 1:
				cur_groups = ALG.collapseGroups(ginfo, hybrid_clade, v);
				gt_groups[gene_num] = cur_groups;
				#print cur_groups;

			outline = "GT-" + str(gene_num) + " to MT-" + str(mul_num) + "\t";
			# Parsing the gene tree.
			if v == -2:
				print 'gt:', gt;

			if not check_nums:
				dup_score, loss_score, maps = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, gt_groups[gene_num], cap, v, check_nums);
			else:
				num_groups, num_fixed = ALG.mulRecon(hybrid_clade, mt, minfo, gt, ginfo, gt_groups[gene_num], cap, v, check_nums);
				outline += str(num_groups) + "\t" + str(num_fixed) + "\t" + str(2**num_groups);
				RC.printWrite(outfilename, v, outline);
				continue;
			# The call of the reconciliation algorithm! On the current gene tree with the current MUL-tree.

			if dup_score == "OVER":
				outline += "Number of groups (" + str(loss_score) + ") over group cap (-p set to " + str(cap) + ") -- Skipping."
				num_skipped += 1;
			else:
				mut_score = dup_score + loss_score;
				outline = outline + str(dup_score) + "\t" + str(loss_score) + "\t" + str(mut_score);
				mul_dict[mul_num][2] += mut_score;
			RC.printWrite(detoutfilename, v, outline);
			# Output stats for the current reconciliation.

		if not check_nums:
			RC.printWrite(detoutfilename, v, "Number of trees skipped: " + str(num_skipped));
			RC.printWrite(detoutfilename, v, "Total parsimony score for MT-" + str(mul_num) + ": " + str(mul_dict[mul_num][2]));
			RC.printWrite(detoutfilename, v, "# ---------------------------");

			RC.printWrite(outfilename, 0, "MT-" + str(mul_num) + "\t" + copy_node + "\t" + RT.mulPrint(mt, hybrid_clade) + "\t" + str(mul_dict[mul_num][2]));

		# Output total score for the current MUL-tree

		mul_num = mul_num + 1;
		group_flag = 0;

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

		RC.printWrite(detoutfilename, 0, "The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + min_tree);
		RC.printWrite(outfilename, 0, "# ---------")		
		RC.printWrite(detoutfilename, 0, "Score = " + str(min_score));
		

RC.printWrite(outfilename, 0, "# The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" + min_tree);
RC.printWrite(outfilename, 0, "# Score = " + str(min_score));
RC.printWrite(detoutfilename, 1, "\nDone!\n");
endtime = time.time();
totaltime = endtime - starttime;
RC.printWrite(outfilename, 1, "# Total execution time: " + str(round(totaltime,3)) + " seconds.");
RC.printWrite(outfilename, 1, "# =========================================================================");





