import sys 
import os 
import multiprocessing as mp
import grampa_lib.reconcore as RC

#############################################################################

def optParse(globs):
# This function handles the command line options and prepares the output directory and files.
# Defaults are set in params.py
    try:
        import psutil
        globs['psutil'] = True;
    except:
        globs['psutil'] = False;
    # Check if psutil is installed for memory usage stats.

    try:
        import argparse;
    except:
        print("\n*** ERROR: Your installation of Python is missing the argparse module. Please try a different version of Python (2.7 or later) or install the module.\n")
        sys.exit();
    # First check if the argparse module is installed. If not, the input options cannot be parsed.

    parser = argparse.ArgumentParser(description="GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis.");

    parser.add_argument("-s", dest="spec_tree", help="A file or string containing a ROOTED, bifurcating, newick formatted species tree in newick format on which to search for polyploid events.");
    parser.add_argument("-g", dest="gene_input", help="A file containing one or more ROOTED, bifurcating (no polytomies), newick formatted gene trees to reconcile. The labels in the gene tree MUST end with '_[species name]' and contain no other underscores. Also accepts a single tree string.");
    parser.add_argument("-h1", dest="h1_spec", help="A space separated list of species labels or internal nodes that define the polyploid clade. Example: 'x,y,z y,z' or '2 4'", default=False);
    parser.add_argument("-h2", dest="h2_spec", help="A space separated list of species labels or internal node labels that make up the clade that you wish to place the second polyploid clade sister to. Example: 'c'", default=False);
    parser.add_argument("-c", dest="group_cap", help="The maxmimum number of groups to consider for any gene tree. Default: 8. Max value: 18.", type=int, default=8);
    parser.add_argument("-o", dest="output_dir", help="Output directory name. Default: grampa-[current date]-[current time]", default="grampa-" + RC.getLogTime());
    parser.add_argument("-f", dest="output_prefix", help="A prefix to add to the beginning of all output files created. Default: grampa", default="grampa");
    parser.add_argument("-p", dest="processes", help="The number of processes GRAMPA should use. Default: 1.", type=int, default=1);
    parser.add_argument("-v", dest="verbosity", help="An option to control the amount of output printed to the screen. 0: print nothing. 1: print only some info at the start. 2: print all log info to screen. 3 (default): print all log info to the screen as well as progress updates for certain steps.", type=int, default=3);
    parser.add_argument("-r", dest="groups_dir", help=argparse.SUPPRESS );
    # help="If you've done an identical run of GRAMPA before, you can specify the path to that runs 'groups_dir' here to avoid calculating the groups again."
    # User params

    parser.add_argument("--multree", dest="mul_tree_input", help="Set this option if your input species tree is a MUL-tree", action="store_true");
    parser.add_argument("--labeltree", dest="label_opt", help="If this flag is set, the program will read your species tree and simply print it out with the internal nodes labeled.", action="store_true");
    parser.add_argument("--buildmultrees", dest="mul_opt", help="Use this along with -s and possibly -h1 and -h2 to simply build MUL-trees from those options.", action="store_true");
    parser.add_argument("--numtrees", dest="num_mul_opt", help="Use this along with -s and possibly -h1 and -h2 to simply count the number of MUL-trees from those options.", action="store_true");
    parser.add_argument("--checknums", dest="check_nums", help="Use this flag in conjunction with all other options to check the number of nodes, groups, and combinations for each gene tree and MUL-tree. In general, gene trees with more than 15 groups to map take a very long time to reconcile.", action="store_true");
    parser.add_argument("--no-st", dest="no_st_opt", help="Set this to only perform reconciliations to the input MUL-trees and NOT the singly-labeled tree.", action="store_true");
    parser.add_argument("--st-only", dest="only_st_opt", help="Set this to only perform reconciliations to the singly-labeled tree.", action="store_true");
    parser.add_argument("--maps", dest="maps_opt", help="Output the maps for each reconciliation in the detailed output file.", action="store_true");
    
    parser.add_argument("--info", dest="info_flag", help="Print some meta information about the program and exit. No other options required.", action="store_true", default=False);
    parser.add_argument("--overwrite", dest="ow_flag", help="Set this to overwrite existing files.", action="store_true", default=False);
    parser.add_argument("--appendlog", dest="append_log_flag", help="Set this to keep the old log file even if --overwrite is specified. New log information will instead be appended to the previous log file.", action="store_true", default=False);
    # User options
    
    parser.add_argument("--tests", dest="test_opt", help="Use 'grampa.py --tests' the first time you run grampa to run through all the options with pre-set input files.", action="store_true");
    # Test script option

    #parser.add_argument("--orthologies", dest="orth_opt", help=argparse.SUPPRESS, action="store_true");
    # help="BETA OPTION: When set, GRAMPA will try to disecern the relationships of polyploid genes in each gene tree (ie paralog vs homoeolog). This method is still under development and may not yet return reliable results!"
    # Beta options

    parser.add_argument("--norun", dest="norun", help=argparse.SUPPRESS, action="store_true", default=False);
    parser.add_argument("--stats", dest="stats_opt", help=argparse.SUPPRESS, action="store_true");
    parser.add_argument("--simpson", dest="s_opt", help=argparse.SUPPRESS, action="store_true");
    # Run mode options

    args = parser.parse_args();
    # The input options and help messages

    warnings = [];
    # List of warnings to print after logfile is created

    globs['call'] = " ".join(sys.argv);
    # Save the program call for later

    ####################

    if args.s_opt:
        RC.simpson();
        sys.exit();
    # ???

    if args.info_flag:
        globs['info'] = True;
        globs['log-v'] = -1;
        startProg(globs);
        return globs;
    # Parse the --info option and call startProg early if set

    if args.norun:
        globs['norun'] = True;
        globs['log-v'] = -1;
    # Parse the --norun option to just parse input info and exit
    
    if args.test_opt:
        RC.testPrep();
        sys.exit();
    # Call of the tests script if --tests is set.

    ## Check run mode options.
    ####################

    if args.spec_tree == None:
        RC.errorOut("OP1", "-s must be specified.", globs);
    globs['spec-tree-input'] = args.spec_tree;
    if not os.path.isfile(globs['spec-tree-input']):
        globs['st-input-type'] = "str";
        globs['st-str'] = globs['spec-tree-input'];
    # The one option required for every single run of GRAMPA is the species tree. This checks that something has been input.
    # Don't check the filename here because a tree string could be input as well

    if args.h1_spec:
        globs['h1-input'] = args.h1_spec;
    if args.h2_spec:
        globs['h2-input'] = args.h2_spec;

    ## Species tree options
    ####################

    if args.check_nums:
        globs['check-nums'] = True;

    if args.maps_opt:
        globs['maps-opt'] = True;

    ####################

    if args.label_opt:
        print("\n*** Message: --labeltree is set to True! Just labeling your species (-s) tree. All other options will be ignored!\n");
        import grampa_lib.spec_tree as ST
        globs = ST.readSpecTree(globs, label_opt=True);
    if args.mul_opt:
        globs['mul-opt'] = True;
        print("\n*** Message: --buildmultrees is set to True! Just printing out all the MUL-trees you requested. All options but -s, -h1, and -h2 will be ignored!\n");
        #globs.cap = 0;
    if args.num_mul_opt:
        globs['num-opt'] = True;
        print("\n*** Message: --numtrees is set to True! Just printing out the number of MUL-trees you requested. All options but -s, -h1, and -h2 will be ignored!\n");
        return(globs);
        #globs.cap = 0;
    # Checking if the --labeltree or --buildmultrees options have been set.

    ####################

    if not args.mul_opt and not args.num_mul_opt and args.gene_input == None:
        RC.errorOut("OP2", "-g must be specified.", globs);
    globs['gt-input'] = args.gene_input;
    if not os.path.isfile(globs['gt-input']):
        globs['gt-input-type'] = "str";
        globs['gt-strs'][1] = globs['gt-input'];
    # The next important input for doing reconciliations is the gene tree file. This checks that something has been input.

    ####################

    if args.no_st_opt and args.only_st_opt:
        RC.errorOut("OP3", "Only one of --no-st and --st-only can be specified.", globs);
    # Check if both --no-st and --st-only have been specified and error if so
    
    if (args.no_st_opt or args.only_st_opt) and args.mul_tree_input:
        warnings.append("# WARNING: With a MUL-tree as input (--multree), --no-st and --st-only are ignored!");
    # Warn if input spec type is MUL tree, in which case these options are ignored
    
    if args.no_st_opt:
        globs['lca-opt'] = "no-st";
    elif args.only_st_opt:
        globs['lca-opt'] = "st-only";
        globs['cap'], args.check_nums = 0, False;
    # Check if only one has been specified

    ## Parse --no-st and --st-only options
    ###############

    if args.mul_tree_input:
        globs['mul-input-flag'] = True;
        if (args.h1_spec != False or args.h2_spec != False):
            warnings.append("# WARNING: With a MUL-tree as the input species tree (--multree) input for -h1 and -h2 are not required and will be ignored.");
    # This checks the input species tree type.

    ###############

    if not any([args.only_st_opt, args.mul_opt, args.num_mul_opt, args.label_opt]):
        globs['cap'] = RC.isPosInt(args.group_cap);
        if not globs['cap']:
            RC.errorOut("OP4", "-c must be a positive integer.", globs);
        elif globs['cap'] > 18:
            RC.errorOut("OP4", "-c cannot be set higher than 18.", globs);
        elif globs['cap'] >= 15:
            warnings.append( "# Warning! With -c set to 15 or higher, some gene trees may take a very long time to reconcile!\n");
    # Checking the group cap.

    ###############

    if args.verbosity not in [-1,0,1,2,3]:
        globs['log-v'] = 3;
        RC.errorOut("OP4", "-v must take values of 0, 1, 2, or 3", globs);
    elif not globs['norun']:
        globs['log-v'] = args.verbosity;
        if globs['log-v'] < 3:
            globs['full-updates'] = False;
    # Checking the verbosity option.

    ###############

    globs['num-procs'] = RC.isPosInt(args.processes);
    if not globs['num-procs']:
        RC.errorOut("OP4", "-n must be a positive integer.", globs);

    ###############

    ### Begin output prep block.
    outdir = os.path.normpath(args.output_dir);
    # Initialize output directory and files.

    ####################

    if not args.label_opt and not args.num_mul_opt:
        if not args.output_dir:
            globs['outdir'] = "grampa-out-" + globs['startdatetime'];
        else:
            globs['outdir'] = args.output_dir;

        if not args.ow_flag and os.path.exists(globs['outdir']):
            RC.errorOut("OP8", "Output directory already exists: " + globs['outdir'] + ". Specify new directory name OR set --overwrite to overwrite all files in that directory.", globs);

        if not os.path.isdir(globs['outdir']) and not globs['norun'] and not globs['info']:
            os.makedirs(globs['outdir']);
        # Main output dir

        if args.output_prefix:
            globs['prefix'] = args.output_prefix;

        globs['output-file'] = os.path.join(globs['outdir'], globs['prefix'] + "-scores.txt");

        ####################

        globs['logfilename'] = os.path.join(globs['outdir'], globs['prefix'] + ".log");

        if not args.append_log_flag and not globs['norun']:
            logfile = open(globs['logfilename'], "w");
            logfile.write("");
            logfile.close();
        # Prep the logfile to be overwritten if --appendlog isn't specified

        if warnings:
            for warning in warnings:
                RC.printWrite(globs['logfilename'], globs['log-v'], warning);
                globs['warnings'] += 1; 
            RC.printWrite(globs['logfilename'], globs['log-v'], "#");
        # Print any warnings here if there were any before the logfile was created  

        ####################

        globs['gt-filtered-file'] = os.path.join(globs['outdir'], globs['prefix'] + "-trees-filtered.txt");

        if not globs['lca-opt'] == "st-only":
            if not args.groups_dir:
                globs['pickle-dir'] = os.path.join(globs['outdir'], "groups-dir.pickle");
                if not globs['mul-opt'] and not globs['num-opt']:
                    os.makedirs(globs['pickle-dir'], exist_ok=True);
            elif os.path.isdir(args.groups_dir):
                globs['pickle-dir'] = args.groups_dir;
            else:
                RC.errorOut(6, "Cannot find specified groups directory! (-r)", globs);

        if not globs['mul-opt']:
            globs['checknums-outfile'] = os.path.join(globs['outdir'], globs['prefix'] + "-checknums.txt");
            # Filtered gene trees file, groups file, and checknums file.

            if not globs['check-nums']:
                globs['detailed-outfile'] = os.path.join(globs['outdir'], globs['prefix'] + "-detailed.txt");
                globs['dup-count-outfile'] = os.path.join(globs['outdir'], globs['prefix'] + "-dup-counts.txt");               
            # If --checknum is not set, we have to prepare detailed and duplication count output files.

                #if globs.orth_opt:
                #    globs.labeled_tree_file = os.path.join(outdir, args.output_prefix + "_labeled_trees.txt");
                #    globs.orth_file_name = os.path.join(outdir, args.output_prefix + "_orthologies.txt");
                # The orthology output files.

    ####################

    startProg(globs);
    # After all the essential options have been set, call the welcome function.
    
    globs['pool'] = mp.Pool(processes=globs['num-procs']);
    if globs['psutil']:
        globs['pids'] = [psutil.Process(os.getpid())];
    # Get the starting process ids to calculate memory usage throughout.

    return globs;

#############################################################################

def startProg(globs):
    # A nice way to start the program.
    if globs['log-v'] != 0:
        start_v = 3;
        print("#");
    else:
        start_v = 1;

    RC.printWrite(globs['logfilename'], start_v, "# =========================================================================");
    RC.printWrite(globs['logfilename'], start_v, "# Welcome to GRAMPA -- Gene tree reconciliations with MUL-trees.");
    RC.printWrite(globs['logfilename'], start_v, "# Version " + globs['version'] + " released on " + globs['releasedate']);
    RC.printWrite(globs['logfilename'], start_v, "# GRAMPA was developed by " + globs['authors']);
    RC.printWrite(globs['logfilename'], start_v, "# Citation:      " + globs['doi']);
    RC.printWrite(globs['logfilename'], start_v, "# Website:       " + globs['http']);
    RC.printWrite(globs['logfilename'], start_v, "# Report issues: " + globs['github']);
    RC.printWrite(globs['logfilename'], start_v, "#");
    RC.printWrite(globs['logfilename'], start_v, "# The date and time at the start is:  " + RC.getDateTime());
    RC.printWrite(globs['logfilename'], start_v, "# Using Python executable located at: " + globs['pyexe']);
    RC.printWrite(globs['logfilename'], start_v, "# Using Python version:               " + globs['pyver'] + "\n#");    
    RC.printWrite(globs['logfilename'], start_v, "# The program was called as:          " + globs['call'] + "\n#");

    if globs['info']:
        return;
    # If --info is set, return after printing program info

    #######################

    pad = 40;
    opt_pad = 30;
    RC.printWrite(globs['logfilename'], start_v, "# " + "-" * 125);
    RC.printWrite(globs['logfilename'], start_v, "# INPUT/OUTPUT INFO:");

    if globs['st-input-type'] == "file":
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Species tree file:", pad) + globs['spec-tree-input']);
    else:
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Input species tree string:", pad) + globs['spec-tree-input']);

    if not globs['mul-opt']:
        if globs['gt-input-type'] == "file":
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Gene tree file:", pad) + globs['gt-input']);
        else:
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Input gene tree string:", pad) + globs['gt-input']);    

    RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Output directory:", pad) + globs['outdir']);

    if not globs['mul-opt']:
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Score file:", pad) + globs['output-file']);
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Filtered gene trees:", pad) + globs['gt-filtered-file']);
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Check nums file:", pad) + globs['checknums-outfile']);
        #RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Groups directory:", pad) + globs['pickle-dir']);
        if not globs['check-nums']:
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Detailed mapping file:", pad) + globs['detailed-outfile']);
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Duplication count file:", pad) + globs['dup-count-outfile']);
        # Input/Output
        #######################

        RC.printWrite(globs['logfilename'], start_v, "# " + "-" * 125);
        RC.printWrite(globs['logfilename'], start_v, "# OPTIONS INFO:");
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# Option", pad) + RC.spacedOut("Current setting", opt_pad) + "Current action");

        opt_str = "All";
        if globs['h1-input']:
            opt_str = globs['h1-input'];
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# -h1", pad) +
                    RC.spacedOut(opt_str, opt_pad) +
                    "GRAMPA will search these H1 nodes. If none are specified, all nodes will be searched as H1 nodes.");
        # -h1

        opt_str = "All";
        if globs['h2-input']:
            opt_str = globs['h2-input'];
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# -h2", pad) +
                    RC.spacedOut(opt_str, opt_pad) +
                    "GRAMPA will search these H2 nodes. If none are specified, all nodes will be searched as H2 nodes.");        
        # -h2

        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# -c", pad) +
                    RC.spacedOut(str(globs['cap']), opt_pad) +
                    "Gene trees with more than this number of groups/clades with polyploid species for a given h1/h2 combination will be skipped.");
        # -c

        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# -f", pad) +
                    RC.spacedOut(globs['prefix'], opt_pad) +
                    "All output files generated will have this string preprended to them.");
        # -f

        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# -p", pad) +
                    RC.spacedOut(str(globs['num-procs']), opt_pad) +
                    "GRAMPA will use this number of processes for LCA mapping.");
        # -p

        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# -v", pad) +
                    RC.spacedOut(str(globs['log-v']), opt_pad) +
                    "Controls the amount of info printed to the screen as GRAMPA is running.");
        # -v

        opt_str = "The tree input with -s will be read as singly-labeled tree.";
        if globs['mul-input-flag']:
            opt_str = "The tree input with -s will be read as a MUL-tree.";
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# --multree", pad) +
                    RC.spacedOut(str(globs['mul-input-flag']), opt_pad) + opt_str);
        # --multree

        opt_str = "GRAMPA will count groups to filter gene trees and then perform reconciliations.";
        if globs['check-nums']:
            opt_str = "GRAMPA will count groups to filter gene trees and exit.";
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# --checknums", pad) +
                    RC.spacedOut(str(globs['check-nums']), opt_pad) + opt_str);
        # --checknums

        opt_str = "GRAMPA will perform reconciliations to all MUL-trees specified by -h1 and -h2 and the input species tree.";
        if globs['lca-opt'] == "no-st":
            opt_str = "GRAMPA will perform reconciliations to only the MUL-trees specified by -h1 and -h2.";
        if globs['lca-opt'] == "st-only":
            opt_str = "GRAMPA will perform reconciliations to only the input species tree.";
        RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# --no-st, --st-only", pad) +
                    RC.spacedOut(str(globs['lca-opt']), opt_pad) + opt_str);
        # --no-st and --st-only

        if not globs['check-nums']:
            opt_str = "GRAMPA will only output duplication and loss counts in the detailed output file.";
            if globs['maps-opt']:
                opt_str = "GRAMPA will output node mappings for the lowest scoring tree in the detailed output file.";
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# --maps", pad) +
                        RC.spacedOut(str(globs['maps-opt']), opt_pad) + opt_str);
        # --maps

        if globs['overwrite']:
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# --overwrite", pad) +
                        RC.spacedOut("True", opt_pad) +
                        "GRAMPA will OVERWRITE the existing files in the specified output directory.");
        # --overwrite

        if globs['norun']:
            RC.printWrite(globs['logfilename'], start_v, RC.spacedOut("# --norun", pad) +
                        RC.spacedOut("True", opt_pad) +
                        "ONLY PRINTING RUNTIME INFO.");
            RC.printWrite(globs['logfilename'], start_v, "# " + "-" * 125);
        # Reporting the norun option

    if globs['log-v'] == 1:
        RC.printWrite(globs['logfilename'], start_v, "# " + "-" * 125);
        RC.printWrite(globs['logfilename'], start_v, "# " + RC.getDateTime() + " INFO: Starting GRAMPA. With -v 1 set, no more information will be printed to the screen until the end of the run.");
        RC.printWrite(globs['logfilename'], start_v, "# " + RC.getDateTime() + " INFO: Check log file for progress updates: " + globs['logfilename']);

    # Other options
    #######################