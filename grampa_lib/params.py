#############################################################################
# This file holds some global variables for some of the input options.
# These global parameters should be read only -- they are not modified anywhere
# else in the code except when reading the input options.
#
# This dictionary will also be used to track output.
#############################################################################

import sys
import os
import timeit
import grampa_lib.reconcore as RC

#############################################################################

class StrictDict(dict):
# This prevents additional keys from being added to the global params dict in
# any other part of the code, just to help me limit it's scope
# https://stackoverflow.com/questions/32258706/how-to-prevent-key-creation-through-dkey-val
    def __setitem__(self, key, value):
        if key not in self:
            raise KeyError("{} is not a legal key of this StrictDict".format(repr(key)));
        dict.__setitem__(self, key, value);

#############################################################################

def init():
    globs_init = {
        'version' : '1.4.0',
        'releasedate' : "March 2023",
        'authors' : "Gregg Thomas, S. Hussain Ather, Matthew Hahn",
        'doi' : 'https://doi.org/10.1093/sysbio/syx044',
        'http' : 'https://gwct.github.io/grampa/',
        'github' : 'https://github.com/gwct/grampa',
        'starttime' : timeit.default_timer(),
        'startdatetime' : RC.getOutTime(),
        # Meta info

        'pyver' :  ".".join(map(str, sys.version_info[:3])),
        'pyexe' : os.path.realpath(sys.executable),
        # System info

        'call' : "",
        # Script call info

        'spec-tree-input' : False,
        'mul-input-flag' : False,
        'st-input-type' : "file",
        'st-str' : False,
        'st' : False,
        'tips' : False,
        'parsed-st-str' : False,
        'h1-input' : False,
        'h2-input' : False,

        'h1-clades' : False,
        'h1-nodes' : False,
        'h2-clades' : False,
        'h2-nodes' : False,

        'mul-trees' : {},
        # mul_trees -> mul_num : [MUL-tree string, MUL-tree dict, hybrid clade species, hybrid node, copy node,
        #  an empty list to add gene tree groups to later]
        
        'min-num' : False,
        'min-score' : False,
        'min-tree' : False,

        'gt-input' : False,
        'gt-input-type' : "file",
        'gt-strs' : {},
        'gt' : {},
        'gt-filtered' : {},
        # gene_trees_filtered will be a list of lists. One list for each gene tree. If the gene tree passes all
        # filters, it will be [gene tree string, gene tree dict]. Otherwise it will be [Filter message].
        # After filtering is done, the filtered trees are removed from the dict
        
        'cap' : 8,
        'lca-opt' : "default",
        'mul-opt' : False,
        'num-opt' : False,
        'check-nums' : False,
        'maps-opt' : False,
        'orth-opt' : False,

        'labeled-tree-file': False,
        'orth-file-name' : False,

        'stats' : False,

        'outdir' : False,
        'prefix' : "grampa",
        'output-file' : False,
        'output-headers' : ['mul.tree', 'h1.node', 'h2.node', 'score', 'labeled.tree'],

        'overwrite' : False,
        
        'detailed-outfile' : False,
        'detailed-headers' : ['mul.tree', 'gene.tree', 'dups', 'losses', 'total.score', 'maps'],

        'dup-count-outfile' : False,
        'dup-count-headers' : ['mul.tree', 'node', 'dups'],

        'checknums-outfile' : False,
        'checknums-headers' : ['mul.tree', 'gene.tree', 'groups', 'fixed', 'combinations', 'over.cap.filtered'],
        'gt-filtered-file' : False,

        'pickle-dir' : False,

        'num-procs' : 1,
        'pool' : False,
        # Number of processes to use and their pool

        'info' : False,
        'quiet' : False,
        # Other user options

        'warnings' : 0,
        'pad' : 82,
        'endprog' : False,
        'exit-code' : 0,
        'log-v' : 1,
        'full-updates' : True,
        'progstarttime' : 0,
        'stepstarttime' : 0,
        'pids' : "",
        'psutil' : False,
        'info' : False,
        'norun' : False,
        'debug' : False,
        'nolog' : False,
        # Internal stuff
    }

    globs_init['logfilename'] = "grampa-" + globs_init['startdatetime'] + ".errlog";
    # Add the runtime to the error log file name.

    globs = StrictDict(globs_init);
    # Restrict the dict from having keys added to it after this

    return globs;

#############################################################################