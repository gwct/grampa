# GRAMPA
### Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis

## Author
#### Gregg Thomas

## About

#### GRAMPA is a program to place polyploidy events on a phylogeny and count duplications and losses in the presence of polyploidy.

## Version History
#### This is version 1.1, released March 22, 2017
* Added reconciliation to the input singly-labeled tree. This allows for comparisions of scenarios of no polyploidy (singly-labeled tree) vs. polyploidy (MUL-trees).
* The -d option has been added to control whether to do reconciliation to the singly-labeled tree. See below for usage.
* Now runs the --checknums option for all reconciliation runs in order to filter trees that are over the group cap for any MUL-tree considered. This ensures that the same trees are run for every MUL-tree.
* Added the --maps option to print the maps for each reconciliation.
* The --multree option to simply build MUL-trees with given singly-labeled tree and h1 and h2 nodes has been renamed --buildmultree.
* The -t option has been renamed with --multree.
* The max group cap has been raised to 18.
* -o now specifies an output DIRECTORY instead of a file. GRAMPA will make this directory for you.
* By default, all output files will be placed in the directory specified with -o with "grampa_" as their prefix (i.e. grampa_out.txt, grampa_det.txt). -p has been added so users can change this prefix.
* Major restructure of functions and libraries.

###### Version 1.0 (Summer, 2016): First release and implementation of MUL-tree reconciliation algorithm.

## Citation

The paper is currently available as a pre-print on bioRxiv: http://biorxiv.org/content/early/2016/06/10/058149

Thomas GWC, Ather SH, Hahn MW. (pre-print) Gene-tree reconciliation with MUL-trees to resolve polyploid events.

## Installation

Simply download the program and run it. You may want to add the GRAMPA folder to your $PATH variable for ease of use.
### The only dependency is Python 2.7 or higher

## Usage

The first thing you should do when you try to run GRAMPA is make sure everything is working with some test files. You can do this easily by running the --tests command:

`python grampa.py --tests`

If all tests pass, then you're good to go! Basic usage in a real case would be:

`python grampa.py -s [species tree file] -g [gene trees file] -o [output directory]`

This would perform a full search for the optimal (lowest-scoring) MUL-tree on the input species tree.

### Input

There are two main inputs for the program. 

1. A file containing a Newick formatted **rooted** species tree (`-s`)
2. A file containing a list of **rooted** gene trees in Newick format (one tree per line). (`-g`)

The labels of the gene tree MUST be formatted such that they end with _[species label], where [species_label] corresponds to a tip label in the species tree.

The other two important input options are `-h1` and `-h2`. See below for their description.

### Output

All output files will be placed in the output directory specified by `-o`

For most runs, GRAMPA creates three main output files: 

1. grampa_out.txt: The main output gives the total reconciliation score for each MUL-tree considered, along with information about the trees. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.

2. grampa_checknum.txt: GRAMPA must calculate how many combinations of maps there are for each MUL-tree/gene tree combo and filter out those that are over the group cap in any combo. The number of groups is recorded in this file.

3. grampa_det.txt: The secondary output file is a detailed output describing the reconciliation scores from each gene tree to each MUL-tree. It has the same name as the file specified with `-o`, but with the suffix '_det' added to the end.

### Options

| Option | Description | 
| ------ | ----------- |
| -s | A file containing a bifurcating, rooted species tree in Newick format. This tree can either be singly-labeled or MUL. |
| -g | A file containing one or more rooted, Newick formatted gene trees. |
| -h1 | A space separated list of nodes to search as the polyploid clade. Only used with `-t s`. If nothing is entered all nodes will be considered. |
| -h2 | A space separated list of nodes to search as possible parental lineages for all nodes specified with `-h1`. If nothing is entered all possible nodes for the current `h1` will be considered. |
| -d | An option to specify whether to do reconciliations to MUL-trees only (0), the singly-labeled tree only (1), or both (2). Default: 2 |
| -c | The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 15 |
| -o | Output directory name. If the directory is not present, GRAMPA will created it for you. |
| -p | By default, all output files created by GRAMPA will have the prefix 'grampa_'. You can specify a different prefix with this option. |
| -v | Control the amount of output printed to the screen. Print all output (1) or just a progress bar and some log info (0). Default: 1 |
| --multree | Set this flag if your input species tree is a MUL-tree. |
| --labeltree | The program will simply label your input species tree. |
| --buildmultrees | Build MUL-trees given `-s`, `-h1`, and `-h2`. |
| --checknums | If this flag is entered, the program will just calculate the number of groups per gene tree and exit. No reconciliations will be done. |
| --maps | Output the node maps for each reconciliation in addition to the scores. The maps will be placed in the detailed output file. |
| --tests | Run the tests script |

### Detailed options

#### `-s` : A file containing a rooted, newick formatted species tree. This tree can be either single-labeled or MUL. 

Entering a singly-labeled tree means you wish to search for the most parsimonious polyploidy scenario. GRAMPA will build MUL-trees based on this singly-labeled tree and calculate reconciliation scores. You can specify the range of MUL-trees to build with the `-h1` and `-h2` options.

Example singly-labeled tree: `(((a,(x,(y,z))),b),(c,d))`

Entering a MUL-tree is the equivalent of entering a standard tree and specifying a single H1 and single H2 node. It represents a single scenario of polyploidy and should be used if you wish to count the number of duplications and losses on gene trees given that scenario. **NOTE: If a MUL-tree is entered, the `--multree` flag must be set.**

Example MUL-tree: `((((a,(x,(y,z))),b),(x,(y,z))),(c,d))`

#### `-g` : A file containing newick formatted gene trees.

This file should contain one or more rooted, Newick formatted gene trees, with one tree per line in the file. 

**The labels of the gene trees MUST end with `_[species label]`**

Where `[species label]` matches a tip label in the species tree. This is necessary so that we can initialize the mappings correctly.

#### `-h1` and `-h2` : GRAMPA's search parameters

H1 and H2 are nodes in the species tree that define how to build a MUL-tree. H1 is the node that represents the polyploid clade. The subtree rooted at H1 and the branch that H1 subtends will be copied onto the branch that H2 subtends:

<p align="center"><img src="http://cgi.soic.indiana.edu/~grthomas/grampa_portal/img/h1h2ex.png"></p>

In the above example, H1 is node 2 and H2 is node 4 in the species tree. This leads to the MUL-tree on the right.

H1 and H2 can be input in 2 different ways:

`-h1 "2" -h2 "4"` and `-h1 "x,y,z" -h2 "a,b,x,y,z"` are equivalent ways of specifying H1 and H2 in the image above. 

The first way relies on internal node labels. To label your species tree, use the `--labeltree` option (described below). **IMPORTANT: For now, only use node labels as specified by `--labeltree`. Custom labels will not work.**

The second way uses the species that define that node. **Species must be comma delimited.**

**H2 cannot be located below H1 in the species tree!** If this occurs, GRAMPA will just tell you that it's not possible and move on.
 
Multiple H1 and H2 nodes can be entered as a space delimited list:

`-h1 "2 4" -h2 "5 6"` and `-h1 "x,y,z a,b,x,y,z" -h2 "c,d a,b,c,d,x,y,z"` are equivalent. 
Entering this means that GRAMPA will first set H1 as node 2 and try both nodes 5 and 6 as H2. Then H1 will be set to node 4 and will try nodes 5 and 6 as H2. 

**If `-h1` and `-h2` are not specified, GRAMPA will search try all possible node combinations of H1 and H2!**

#### `-d` : Reconciliation type option

GRAMPA can do reconciliations to singly-labeled and MUL-trees. If you know a polyploidy event has taken place, you may wish to only reconcile to MUL-trees. However, if you are trying to identify a new polyploidy event, the scores of all MUL-trees considered must be compared to the score of the singly-labeled tree, which represents a scenario of no polyploidy.

To reconcile to MUL-trees only, set `-d 0`

To reconcile to the singly-labeled tree only, set `-d 1`

To reconcile to both the singly-labeled and MUL-trees set `-d 2` (Default)

#### `-c` : The group cap

GRAMPA uses the standard LCA reconciliation algorithm on MUL-trees, meaning that some genes have more than one possible mapping. We get around this by trying ALL possible initial mappings and picking the one with the lowest score. This works, but also means our program has an exponential runtime based on the number of genes from polyploid species in any given gene tree. We get around this in several ways by collapsing and fixing groups (see paper), but there can still be lots of groups. This parameter sets the maximum number of groups to consider for any gene tree. If a gene tree has more than this number of groups, it will be skipped. 

Default is 8 groups, with a max setting of 18.

#### `-o` : Output directory

Grampa creates many output files, so it is easiest just to place them all in a single directory. That directory can be specified with this option, and will be created for you if it doesn't exist. If this option is not specified, the default output directory is "grampa_[date]-[time]".

There are two output files created by GRAMPA. 

Depending on the flags set, GRAMPA creates three main output files. They are described above in the **Output** section.

#### `-p` : Output file prefix

By default, all output files created by GRAMPA will have the prefix 'grampa_'. You can specify a different prefix with this option.
For example, a run with `-p test` will generate the following output files, all within the output directory:

`test_out.txt, test_det.txt, test_checknums.txt`

#### `--multree` : Input tree flag

GRAMPA can accept both singly-labeled and MUL-trees as input. If you input species tree (`-s`) is a MUL-tree, you must set this flag.
A MUL-tree represents a single possible polyploid scenario and it is equivalent to entering a singly-labeled tree with a single H1 and H2 node specified.

#### `--labeltree` : Tree labeling

This option can be used in conjunction with `-s` to simply add internal node labels to a species tree. 

`python grampa.py -s species.tree --labeltree`

If the input in `species.tree` is : `(((a,(x,(y,z))),b),(c,d))`
Then the output tree will be: `(((a,(x,(y,z)<1>)<2>)<3>,b)<4>,(c,d)<5>)<6>`

#### `--buildmultrees` : Build MUL-trees

This option can be used with `-s`, `-h1`, and `-h2` to build MUL-trees from a standard species tree.

`python grampa.py -s species.tree -h1 "2" -h2 "4" -o multree.txt --multree`

If the input in `species.tree` is : (((a,(x,(y,z))),b),(c,d))
Then the output MUL-tree will be: ((((a,(x+,(y+,z+)<1>)<2>)<3>,b)<4>,(x*,(y*,z*)<5>)<6>)<7>,(c,d)<8>)<9>

#### `--checknums` : Group counting

With this set, the program will run normally with the specified options, except no reconciliations will be done. Instead, only the checknums output file will be created and will contain the number of polyploid groups for each gene tree. Use this to decide the best setting for `-c`.

#### `--maps` : Output node mappings

Set this option to output the LCA node mappings along with the reconciliation scores to the detailed (_det.txt) output file.

#### `--tests` : Tests script

This option runs the tests.py script to ensure that all the other options in the program are working.

## License

This file is part of GRAMPA.

GRAMPA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GRAMPA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 






