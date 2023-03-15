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

import sys 
import os
import timeit
import multiprocessing as mp
from functools import partial
import pickle

import grampa_lib.reconcore as CORE
import grampa_lib.mul_recon as ALG
import grampa_lib.opt_parse as OP
import grampa_lib.mul_tree as MT
import grampa_lib.spec_tree as ST
import grampa_lib.gene_tree as GT
import grampa_lib.params as params
import grampa_lib.mul_out as OUT

#############################################################################

def grampa(globs):
    ###########################
    ### Prep

    step_start_time = CORE.report_step(globs, "", "", "", start=True);
    # Initialize the step headers

    ###########################
    ### Species tree

    step = "Reading species tree";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    globs = ST.readSpecTree(globs);
    hybrid_nodes, copy_nodes = "", "";
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: species tree read");
    
    ###########################
    ### Hybrid clades

    if globs['lca-opt'] != "st-only":
        step = "Parsing hybrid clades";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        globs = ST.hInParse(globs);
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: got H nodes");
        
        step = "Counting MUL-trees to be generated";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        num_mul_trees = MT.countMULTrees(globs);
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(num_mul_trees) + " total MUL-trees");
    else:
        num_mul_trees = 0;

    ###########################
    ### Build MUL-trees

    globs['mul-trees'] = { 0 : [globs['st-str'], globs['st'], "", "", ""] };
    # Add the species tree to the MUL-tree dict

    if globs['mul-input-flag']:
         globs['mul-trees'] = { 1 : [globs['st-str'], globs['st'], globs['h1-clades'][0], globs['h1-nodes'][0], globs['h2-nodes'][0], []] };
    # If the input tree is a MUL-tree, we just need to set some variables.
    else:
        step = "Building MUL-trees";
        step_start_time = CORE.report_step(globs, step, False, "In progress...");
        globs = MT.genMULTrees(globs, step_start_time);
        step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['mul-trees'])-1) + " MUL-trees built");
        # genMULTree does all the parsing of the h1 and h2 clades and calls the buildMulTree function.

    ###########################
    ### Gene trees

    step = "Reading gene trees";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");

    skipped_lines = [];

    if globs['gt-input-type'] == "file":
        try:
            gene_num = 1;
            for line in open(globs['gt-input'], "r"):
                globs['gt-strs'][gene_num] = line.strip();
                gene_num += 1;
        except:
            CORE.errorOut("MAIN1", "Error reading gene trees file!", globs);
    # If the gene tree input is a file, try reading the file.

    for result in globs['pool'].imap_unordered(GT.readGeneTree, globs['gt-strs'].items()):
        globs['gt-filtered'][result[0]] = result[1];
        if result[2]:
            skipped_lines.append(result[0]);   
    # Parsing the gene trees to get info about each node in each tree.

    if len(skipped_lines) == len(globs['gt-strs']):
        CORE.errorOut("MAIN2", "Couldn't find any gene trees in your gene tree input file (-g)!", globs);
    # If the input is a file, we assume each line contains one gene tree.

    num_filtered_gt = len(skipped_lines)
    step_start_time = CORE.report_step(globs, step, step_start_time, "Success: " + str(len(globs['gt-strs']) - num_filtered_gt) + " gene trees read");

    if skipped_lines:
        for line_num in skipped_lines:
            CORE.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Line " + str(line_num) + " could not be parsed as a tree and was skipped!");
            globs['warnings'] += 1;

    ###########################
    ### Collapse groups
    if not globs['lca-opt'] == 'st-only':
        step = "Collapsing gene tree groupings";
        step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=globs['full-updates']);
        # Status update

        if globs['mul-input-flag']:
            for mul_num, mul_tree in globs['mul-trees'].items():
                ALG.collapseGroups((mul_num, mul_tree), globs['gt-filtered'], globs['mul-input-flag'], globs['log-v'], globs['pickle-dir'], num_mul_trees);
        ### SERIAL AND MUL-TREE VERSION                
        else:
            for result in globs['pool'].imap_unordered(partial(ALG.collapseGroups, gene_trees_filtered=globs['gt-filtered'], mul_input_flag=globs['mul-input-flag'], v=globs['log-v'], pickle_dir=globs['pickle-dir'], nmt=num_mul_trees), globs['mul-trees'].items()):
                pass;
        ### PARALLEL AND SINGLY LABELED TREE VERSION
        # If the input tree is singly-labeled, do some multi-processing, else if it is a MUL-tree, just run the function.
        ## Feb 14 -- this seems like it could be simplified but I'm not messing with it right now.

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=globs['full-updates']);
    # Status update

    globs, num_gt_over_cap = OUT.checkOut(globs);
    num_filtered_gt += num_gt_over_cap;
    # Filter gene trees over group cap

    gene_trees = OUT.filterOut(globs, num_filtered_gt);
    # Write the filtered trees to a file, or not if no filtering was done.

    if globs['check-nums']:
        CORE.endProg(globs);
    # If --checknums is set, exit the program here.
    CORE.printWrite(globs['logfilename'], globs['log-v'], "# ----------------------------------------");

    ###########################
    ### Reconciliations

    step = "Reconciliation";
    step_start_time = CORE.report_step(globs, step, False, "In progress...", full_update=globs['full-updates']);
    # Status update
    
    all_scores, min_num, min_score, min_maps = {}, '', 9999999, {};
    # Return values for the reconciliations
    # all_scores: dict of scores for all MUL-trees
    # min_num: the number of the lowest scoring MT
    # min_score: the score of the lowest scoring MT
    # min_maps: the detailed mapping of the lowest scoring MT

    for result in globs['pool'].imap_unordered(partial(ALG.mulRecon, gene_trees=globs['gt-filtered'], v=globs['log-v'], pickle_dir=globs['pickle-dir'], nmt=num_mul_trees), globs['mul-trees'].items()):
        all_scores[result[0]] = result[1];
        # Save the score of the current MT in all_scores

        if result[1] < min_score:
            min_num = result[0];
            min_score = result[1];
        # Check if the current score is lower than the lowest so far and update

        del result;
    ## End recon loop

    globs['min-num'] = min_num;
    globs['min-score'] = min_score;
    globs['min-tree'] = globs['mul-trees'][globs['min-num']];

    all_scores = [ (k, v) for k, v in sorted(all_scores.items(), key=lambda item: item[1]) ];
    # Sort the all_scores dict by score and store as a list of tuples

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success", full_update=globs['full-updates']);
    # Status update

    ###########################
    ### Counts for low scoring trees

    step = "Getting maps for lowest scoring MTs";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    lowest_scoring_trees = [ (score_tuple[0], globs['mul-trees'][score_tuple[0]]) for score_tuple in all_scores[:6] ];
    # Get a list of tuples for the 6 lowest scoring MUL-trees
    # [ (mt num, mt string) ]

    lowest_maps = [];
    # A list to store the mappings for each of the lowest scoring MTs

    for low_tree in lowest_scoring_trees:
        lowest_maps.append(ALG.mulRecon(low_tree, globs['gt-filtered'], 0, globs['pickle-dir'], nmt=num_mul_trees, retmap=True))
    ## This gets the mul_tree with the min score along with its maps.

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    ###########################
    ### Output


    step = "Writing detailed output file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    multiple_maps = OUT.detOut(globs, lowest_maps[0]);
    # Gets the detailed mapping info for the lowest scoring tree

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    ##########

    step = "Writing main output file";
    step_start_time = CORE.report_step(globs, step, False, "In progress...");
    # Status update

    OUT.mainOut(globs, all_scores, lowest_maps, multiple_maps);

    step_start_time = CORE.report_step(globs, step, step_start_time, "Success");
    # Status update

    ###########################
    ### Begin orthology prediction block. *BETA*
    if globs['orth-opt']:
        import lib.orth_label as OL
        step = RC.printStep(step, "# " + RC.getDateTime() + " --> STEP " + str(step) + " " + RC.getLogTime() +  ": *BETA* Mapping polyploid genes.");
        min_tree = mul_trees[min_num][0];
        min_clade = mul_trees[min_num][2];
        OL.orthLabel(gene_trees, min_maps, min_tree, min_clade);

    ###########################
#############################################################################

if __name__ == '__main__':
# Main is necessary for multiprocessing to work on Windows.

    start = timeit.default_timer();
    globs = params.init();
    # Get the global params as a dictionary.

    if "-v 0" not in " ".join(sys.argv):
        print("\n" + " ".join(sys.argv) + "\n");
    
    if any(v in sys.argv for v in ["--version", "-version", "--v"]):
        print("# GRAMPA version " + globs['version'] + " released on " + globs['releasedate'])
        sys.exit(0);
    # The version option to simply print the version and exit.
    # Need to get actual degenotate version for this, and not just the interface version.

    #print("#");
    #print("# " + "=" * 100);
    #print(CORE.welcome());
    #if "-h" not in sys.argv:
    #    print("             Gene-tree reconciliation algorithm with MUL-trees for polyploid analysis\n");
    # A welcome banner.
    
    globs = OP.optParse(globs);
    # Getting the input parameters from optParse.

    if globs['info']:
        print("# --info SET. EXITING AFTER PRINTING PROGRAM INFO...\n#")
        sys.exit(0);
    if globs['norun']:
        print("# --norun SET. EXITING AFTER PRINTING OPTIONS INFO...\n#")
        sys.exit(0);
    # Early exit options

    grampa(globs);
    if globs['log-v'] > 0:
        print("# " + CORE.getDateTime() + " INFO: GRAMPA has finished!");
    CORE.endProg(globs);

#############################################################################