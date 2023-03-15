import sys
import re
import grampa_lib.recontree as RT
import grampa_lib.reconcore as RC

#############################################################################
def countMULTrees(globs):
# Function to count the number of possible MUL-trees given a set of h1 and h2 nodes.
    nt = len(globs['st']);
    num_mul_trees = 0;
    for n1 in globs['h1-nodes']:
        # print(n1);
        if globs['st'][n1][2] == 'tip':
            n1_clade = [];
        else:
            n1_clade = RT.getCladeNode(n1, globs['st']);

        ni = 0;
        for n2 in globs['h2-nodes']:
            if n2 not in n1_clade:
                ni += 1;
        num_mul_trees += (ni);

    if globs['num-opt']:
        print();
        print("Total nodes in species tree: " + str(nt));
        print("Total tips in species tree.: " + str(len(globs['tips'])));
        print("H1 nodes...................: " + ",".join(globs['h1-nodes']));
        print("H2 nodes...................: " + ",".join(globs['h2-nodes']));
        print("Possible MUL-trees.........: " + str(num_mul_trees));
        print();
        sys.exit(0)

    return num_mul_trees;

#############################################################################

def genMULTrees(globs, step_start_time):
# genMULTree generates all MUL-trees given a species tree and a set of h1 and h2 nodes.

    if globs['mul-opt']:
        step_start_time = RC.report_step(globs, "Building MUL-trees", step_start_time, "--buildmultrees: exiting after");
        RC.printWrite(globs['logfilename'], globs['log-v'], "\t".join(globs['output-headers'][:-1]));
    # Status update update if --buildmultrees set

    mul_num = 1;

    for h1_node in globs['h1-nodes']:
        h1_clade = set(RT.getClade(h1_node, globs['st']));

        for h2_node in globs['h2-nodes']:
            mt_unlabel = buildMultree(h1_node, h2_node, globs['parsed-st-str'], globs['st']);        
            # Building the MUL-tree by passing the species tree to the buildMultree function
            # Input is one node at which to copy the subtree (hybrid_node)
            # and one node at which to place the copy (copy_node)

            if mt_unlabel == "NULL":
                #print("\n*** Warning: H2 node (" + copy_node + ") within hybrid subtree rooted at H1 node (" + hybrid_node + ") or at root of tree. Not building MUL-tree for this combination.\n");
                #print("# ---------------------------");
                continue;
            # Nodes within the hybrid subtree cannot be copy nodes.

            minfo, mt = RT.treeParse(mt_unlabel);
            # Reading the MUL-tree as usual with treeParse.

            globs['mul-trees'][mul_num] = [mt, minfo, h1_clade, h1_node, h2_node]
            #mul_trees[mul_num] = [mt, minfo, hybrid_clade, hybrid_node, copy_node];
            
            if globs['mul-opt']:
                outline = str(mul_num) + "\t" + h1_node + "\t" + h2_node + "\t" + mulPrint(mt, h1_clade);
                RC.printWrite(globs['logfilename'], globs['log-v'], outline);
            # Print the tree if --buildmultrees is set
            mul_num += 1;

    # This block builds the MUL-trees and prepares the main mul_dict:
    # mul_trees -> mul_num : [MUL-tree string, MUL-tree dict, hybrid clade species, hybrid node, copy node, an empty list to add gene tree groups to later]
    # mul_results -> mul_num : [an initial cumulative score of 0, { gene_num -> [[score for gt/mt combo, # dups for gt/mt combo, # losses for gt/mt combo, maps, node_dups, node_loss],[ditto if multiple maps with same score]]}

    if globs['mul-opt']:
        RC.endProg(globs);
    ## If --buildmultrees is set, this just prints out the MUL-trees and exits.

    return globs;

#############################################################################

def buildMultree(h, p, tree, tree_info):
# This function builds a MUL-tree from a normally labeled species tree.
# Input:
#        h - a node that will have its subtree copied into the MUL-tree
#        p - a node at which the copied subtree will be placed
#        tree - the original species tree WITH internal nodes labeled
#        tree_info - the tree info dictionary (returned by treeparse) from the original tree

    if tree_info[h][2] == 'tip':
        hybrid = h;
    else:
        hybrid = RT.getSubtree(h,tree);
        hybrid = re.sub('<[\d]+>','',hybrid);
    # Gets the subtree of the hybrid node

    if tree_info[p][2] == 'tip':
        copy = p;
    else:
        copy = RT.getSubtree(p,tree);
        copy = re.sub('<[\d]+>','',copy);
    # Gets the subtree of the copy node.

    if (copy in hybrid and copy != hybrid):# or tree_info[p][3] == 'root':
        return "NULL";
    # Copy nodes shouldn't be within the hybrid subtree and shouldn't be at the root.

    if tree_info[h][2] != 'tip':
        for node in tree_info:
            if node in hybrid and tree_info[node][2] == 'tip':
                tree = tree.replace(node, node+"_1");
                copy = copy.replace(node, node+"_1")
    # Some re-labeling if necessary.
    tree = re.sub('<[\d]+>','',tree);

    if tree_info[h][2] == 'tip':
            mul_clade = "(" + copy + "," + hybrid + "*)";
    else:
            mul_clade = "(" + copy + "," + hybrid + ")";

    mul_tree = tree.replace(copy,mul_clade);
    # Combines the clades and replaces the copy clade in the original tree to create the MUL-tree.

    mul_tree = re.sub('<[\d]+>','',mul_tree);
    if tree_info[h][2] != 'tip':
        for node in tree_info:
            if node in hybrid and tree_info[node][2] == 'tip':
                mul_tree = mul_tree.replace(node, node+"*");
        mul_tree = mul_tree.replace("*_1", "");
    # Some relabling of the hybrid species.

    return mul_tree;

#############################################################################

def mulPrint(mul_tree, hybrid_clade):
# For a given MUL-tree, relabels the hybrid clade species to distinguish within tree viewers.

    for spec in hybrid_clade:
        mul_tree = re.sub(spec + '(?!=\*)', spec + '+', mul_tree);
        mul_tree = mul_tree.replace("+*", "*");
    return mul_tree;

#############################################################################