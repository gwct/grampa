import sys, os, subprocess, recontree as RT, reconcore as RC

#############################################################################

def optParse(errorflag):
# This function handles the command line options and prepares the output directory and files.
	try:
		import argparse;
	except:
		print("*** ERROR: Your installation of Python is missing the argparse module. Please try a different version of Python (2.7 or later), or install the module.")
		sys.exit();
	# First check if the argparse module is installed. If not, the input options cannot be parsed.

	parser = argparse.ArgumentParser(description="GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.");

	parser.add_argument("-s", dest="spec_tree", help="A file or string containing a ROOTED, bifurcating, newick formatted species tree in newick format on which to search for polyploid events.");
	parser.add_argument("-g", dest="gene_input", help="A file containing one or more ROOTED, bifurcating, newick formatted gene trees to reconcile. The labels in the gene tree MUST end with '_[species name]' and contain no other underscores.");
	parser.add_argument("-h1", dest="h1_spec", help="A space separated list of species labels that make up the polyploid clade. Example: 'x,y,z y,z' Or a space separated list of nodes that make up the polyploid clade. Example '2 4'", default=False);
	parser.add_argument("-h2", dest="h2_spec", help="A space separated list of species labels or internal node labels that make up the clade that you wish to place the second polyploid clade sister to. If spec tree type (-t) is m, this option can be ignored. Example: 'c'", default=False);
	parser.add_argument("-c", dest="group_cap", help="The maxmimum number of groups to consider for any gene tree. Default: 8. Max value: 15.", type=int, default=8);
	parser.add_argument("-d", dest="lca_standard", help="0: Do not do reconciliation to the singly-labeled input tree. 1: ONLY do reconciliation to the singly-labeled input tree. 2 (default): Do reconciliation to the singly-labeled input tree along with to desired MUL-trees.", type=int, default=2);
	parser.add_argument("-o", dest="output_dir", help="Output directory name. Default: grampa-[current date]-[current time]", default="grampa-" + RC.getLogTime());
	parser.add_argument("-p", dest="output_prefix", help="A prefix to add to the beginning of all output files created. Default: grampa", default="grampa");
	# parser.add_argument("-r", dest="in_groups", help="If you have the groups computed from a previous run (groups.txt), input that file here to skip computing them again.", default="");
	parser.add_argument("-v", dest="verbosity", help="An option to control the amount of output printed to the screen. 0: print only a progress bar and some info. 1: print some detailed output for each reconciliation (this detailed output is also available by default in the _det output file). Default: 1", type=int, default=1);
	parser.add_argument("--multree", dest="spec_tree_type", help="Set this option if your input species tree is a MUL-tree", action="store_true");
	parser.add_argument("--labeltree", dest="label_opt", help="If this flag is set, the program will read your species tree and simply print it out with the internal nodes labeled.", action="store_true");
	parser.add_argument("--buildmultrees", dest="mul_opt", help="Use this along with -s and possibly -h1 and -h2 to simply build MUL-trees from those options.", action="store_true");
	parser.add_argument("--checknums", dest="check_num", help="Use this flag in conjunction with all other options to check the number of nodes, groups, and combinations for each gene tree and MUL-tree. In general, gene trees with more than 15 groups to map take a very long time.", action="store_true");
	parser.add_argument("--maps", dest="maps_opt", help="Output the maps for each reconciliation in the detailed output file.", action="store_true");
	parser.add_argument("--orthologies", dest="orth_opt", help=argparse.SUPPRESS, action="store_true");
	# help="BETA OPTION: When set, GRAMPA will try to disecern the relationships of polyploid genes in each gene tree (ie paralog vs homoeolog). This method is still under development and may not yet return reliable results!"
	parser.add_argument("--tests", dest="test_opt", help="Use 'grampa.py --tests' the first time you run grampa to run through all the options with pre-set input files.", action="store_true");
	parser.add_argument("--simpson", dest="s_opt", help=argparse.SUPPRESS, action="store_true");

	args = parser.parse_args();
	# The input options and help messages

	if errorflag == 0:
	## The initial call of the function to parse all the options.
		if args.s_opt:
			RC.simpson();
			sys.exit();
		# ???

		if args.test_opt:
			RC.testPrep();
			sys.exit();
		# Call of the tests script if --tests is set.

		if args.spec_tree == None:
			RC.errorOut(1, "-s must be specified.");
		# The one option required for every single run of GRAMPA is the species tree. This checks that something has been input.

		if args.label_opt and args.verbosity != -1:
			print("\n*** Message: --labeltree is set to True! Just labeling your species (-s) tree. All other options will be ignored!\n");
			args.group_cap = 0;
		if args.mul_opt and args.verbosity != -1:
			print("\n*** Message: --buildmultrees is set to True! Just printing out all the MUL-trees you requested. All options but -s, -h1, and -h2 will be ignored!\n");
			args.group_cap = 0;
		# Checking if the --labeltree or --buildmultrees options have been set.

		if args.gene_input == None and not args.label_opt and not args.mul_opt:
			RC.errorOut(2, "-g must be specified.");
		# The next important input for doing reconciliations is the gene tree file. This checks that something has been input.

		if args.lca_standard not in [0,1,2]:
			RC.errorOut(3, "-d must take values of 0, 1, or 2");
		if args.lca_standard == 1:
			print("*** Message: -d set to 1: Only doing reconciliation to the singly-labeled tree.");
			print("*** All settings except -s, -g, -o, -p, and --maps will be ignored!");
			args.group_cap, args.check_num = 0, False;
		# Checking the values of the lca_standard option: 0 = do not do recon. to ST, 1 = ONLY do recon to ST, 2 = do recon to ST and MT.
		# If -d has been set to 1, some other things need to be nullified.

		if args.spec_tree_type:
			if (args.h1_spec != False or args.h2_spec != False):
				print("*** Message: With a MUL-tree as the input species tree (-t m) input for -h1 and -h2 are not required.");
				print("*** Your input for -h1 and -h2 will be ignored.\n");
			if args.lca_standard in [1,2]:
				print("*** Message: With a MUL-tree as the input species tree (-t m) reconciliations cannot be done to a singly-labeled tree.");
				print("*** Message: Your input for -d will be ignored.\n");
				args.lca_standard = 0;
		# This checks the input species tree type.

		if args.group_cap > 18:
			RC.errorOut(4, "For computational reasons, -p should not be set higher than 18.");
		elif args.group_cap >=10:
			print("*** Warning! With -p set to 10 or higher, some gene trees may take a very long time to reconcile!\n");
		# Checking the group cap.

		if args.verbosity not in [-1,0,1]:
			RC.errorOut(5, "-v must take values of either 0 or 1");		
		# Checking the verbosity option.

		### Begin output prep block.
		outdir = os.path.normpath(args.output_dir);
		outfilename, checkfilename, detoutfilename, labeled_tree_file, gene_file_filtered, orth_file_name = "", "", "", "", "", "";
		# Initialize output directory and files.
	

		spec_type = 'm' if args.spec_tree_type else 's';
		# The rest of the code still uses the old --multree (-t) formatting for spec_type.

		if not args.label_opt:
			if os.path.isdir(outdir):
				outdir_suffix = 1;
				outdir_prefix = outdir;
				while os.path.isdir(outdir):
					outdir = outdir_prefix + "-" + str(outdir_suffix);
					outdir_suffix += 1;
			os.system("mkdir " + outdir);
			# Making the output directory

			outfilename = os.path.join(outdir, args.output_prefix + "_out.txt");
			RC.filePrep(outfilename);
			# Preparing the main output file
			gene_file_filtered = os.path.join(outdir, args.output_prefix + "_trees_filtered.txt");
			if not args.mul_opt:
				checkfilename = os.path.join(outdir, args.output_prefix + "_checknums.txt");
				RC.filePrep(checkfilename, "Trees\t# Groups\t# Fixed\t# Combinations\n");
				# Filtered gene trees file, groups file, and checknums file.

				if not args.check_num:
					detoutfilename = os.path.join(outdir, args.output_prefix + "_det.txt");
					RC.filePrep(detoutfilename);
				# If --checknum is not set, we have to prepare detailed output file.

					if args.orth_opt:
						labeled_tree_file = os.path.join(outdir, args.output_prefix + "_labeled_trees.txt");
						RC.filePrep(labeled_tree_file);

						orth_file_name = os.path.join(outdir, args.output_prefix + "_orthologies.txt");
						RC.filePrep(orth_file_name);
					# The orthology output files.
		### End output prep block.

		return args.spec_tree, spec_type, args.gene_input, args.h1_spec, args.h2_spec, args.group_cap, outfilename, \
		checkfilename, detoutfilename, gene_file_filtered, labeled_tree_file, orth_file_name, args.verbosity, args.lca_standard, args.label_opt, \
		args.check_num, args.mul_opt, args.orth_opt, args.maps_opt;

	elif errorflag == 1:
		parser.print_help();
		return;
	## If the function is called from an error, it just prints the help and exits the program.

