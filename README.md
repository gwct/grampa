# GRAMPA
### Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis

## Author
#### Gregg Thomas

## About

#### GRAMPA is a program to date polyploidy events and count duplications and losses in the presence of polyploidy.

#### This is version 1.0
#### The only dependency is Python 2.7 or higher

## Installation

Simply download the program and run it. You may want to add the GRAMPA folder to your $PATH variable for ease of use.

## Usage

The first thing you should do when you try to run GRAMPA is make sure everything is working with some test files. You can do this easily by running the --tests command:

`python grampa.py --tests`

If all tests pass, then you're good to go! Basic usage in a real case would be:

`python grampa.py -s [species tree file] -g [gene trees file] -o [output file name]`

This would perform a full search for the optimal MUL-tree on the species tree.

### Input

There are two main inputs for the program. 

1. A file containing a Newick formatted species tree (`-s`)
2. A file containing a list of gene trees in Newick format (one tree per line). (`-g`)

The labels of the gene tree MUST be formatted such that they end with _[species label], where [species_label] corresponds to a tip label in the species tree.

The other two important input options are `-h1` and `-h2`. See below for their description.

### Output

GRAMPA creates two main output files, both specified with `-o`: 

1. The main output gives the total reconciliation score for each MUL-tree considered, along with information about the trees. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.

2. The secondary output file is a detailed output describing the reconciliation scores from each gene tree to each MUL-tree. It has the same name as the file specified with `-o`, but with the suffix '_det' added to the end.

### Options

| Option | Description | 
| ------ | ----------- |
| -s | A file containing a bifurcating species tree in Newick format. This tree can either be standard or MUL |
| -t | This specifies the type of tree entered in `-s`. m: MUL-tree, s: standard tree. Default: s |
| -g | A file containing one or more Newick formatted gene trees |
| -h1 | A space separated list of nodes to search as the polyploid clade. Only used with `-t s`. If nothing is entered all nodes will be considered |
| -h2 | A space separated list of nodes to search as possible parental lineages for all nodes specified with `-h1`. If nothing is entered all possible nodes for the current `h1` will be considered |
| -c | The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 15 |
| -o | Output file name |
| --checknum | If this flag is entered, the program will just calculate the number of groups per gene tree. No reconciliations will be done |
| --labeltree | The program will simply label your input species tree |
| --multree | Build MUL-trees given `-s`, `-h1`, and `-h2` |
| --tests | Run the tests script |

### Detailed options

#### `-s` : A file containing a newick formatted species tree. This tree can be either standard or MUL. 

Entering a standard tree means you wish to search for the most parsimonious polyploidy scenario. GRAMPA will build MUL-trees based on this standard tree and calculate reconciliation scores. You can specify the range of MUL-trees to build with the `-h1` and `-h2` options.

Example standard tree: `(((a,(x,(y,z))),b),(c,d))`

Entering a MUL-tree is the equivalent of entering a standard tree and specifying a single H1 and single H2 node. It represents a single scenario of polyploidy and should be used if you wish to count the number of duplications and losses on gene trees given that scenario. **NOTE: If a MUL-tree is entered, you MUST set `-t m`.**

Example MUL-tree: `((((a,(x,(y,z))),b),(x,(y,z))),(c,d))`

#### `-t` : Input species tree type.

`-t s` : The tree specified with `-s` is standard. This is the Default setting.
`-t m` : The tree specified with `-s` is a MUL-tree.

#### `-g` : A file containing newick formatted gene trees.

This file should contain one or more gene trees, with one tree per line in the file. 

**The labels of the gene trees MUST end with `_[species label]`**

Where `[species labe]` matches a tip label in the species tree. This is necessary so that we can initialize the mappings correctly.

#### `-h1` and `-h2` : GRAMPA's search parameters

H1 and H2 are nodes in the species tree that define how to build a MUL-tree. H1 is the node that represents the polyploid clade. The subtree rooted at H1 and the branch that H1 subtends will be copied onto the branch that H2 subtends:

<p align="center"><img src="http://cgi.soic.indiana.edu/~grthomas/grampa/images/h1h2ex.png"></p>

In the above example, H1 is node 2 and H2 is node 4 in the species tree. This leads to the MUL-tree on the right.

H1 and H2 can be input in 2 different ways:

`-h1 "2" -h2 "4"` and `-h1 "x,y,z" -h2 "a,b,x,y,z"` are equivalent ways of specifying H1 and H2 in the image above. 

The first way relies on internal node labels. To label your species tree, use the `--labeltree` option (described below). **IMPORTANT: For now, only use node labels as specified by `--labeltree`. Custom labels will not work**

The second way uses the species that define that node. **Species must be comma delimited.**

**H2 cannot be located below H1 in the species tree!** If this occurs, GRAMPA will just tell you that it's not possible and move on.
 
Multiple H1 and H2 nodes can be entered as a space delimited list:

`-h1 "2 4" -h2 "5 6"` and `-h1 "x,y,z a,b,x,y,z" -h2 "c,d a,b,c,d,x,y,z"` are equivalent. 
Entering this means that GRAMPA will first set H1 as node 2 and try both nodes 5 and 6 as H2. Then H1 will be set to node 4 and will try nodes 5 and 6 as H2. 

**If `-h1` and `-h2` are not specified, GRAMPA will search try all possible node combinations of H1 and H2!**

#### `-c` : The group cap

GRAMPA uses the standard LCA reconciliation algorithm on MUL-trees, meaning that some genes have more than one possible mapping. We get around this by trying ALL possible initial mappings and picking the one with the lowest score. This works, but also means our program has an exponential runtime based on the number of genes from polyploid species in any given gene tree. We get around this in several ways by collapsing and fixing groups (see paper), but there can still be lots of groups. This parameter sets the maximum number of groups to consider for any gene tree. If a gene tree has more than this number of groups, it will be skipped. 

Default is 8 groups, with a max setting of 15.

#### `-o` : Output files

There are two output files created by GRAMPA. 

If `-o test.txt` is specified you will get:

1. `test.txt` : The main output file containing information about the run and reconciliation scores for each MUL-tree considered. At the bottom the most parsimonious MUL-tree will be displayed.

2. `test_det.txt` : The detailed output file containing reconciliation scores for each gene tree to each MUL-tree.

#### `--checknum` : Group counting

With this set, the program will run normally with the specified options, except no reconciliations will be done. Instead, the output file will contain the number of polyploid groups for each gene tree. Use this to decide the best setting for `-c`.

#### `--labeltree` : Tree labeling

This option can be used in conjunction with `-s` to simply add internal node labels to a species tree. 

`python grampa.py -s species.tree --labeltree`

If the input in `species.tree` is : `(((a,(x,(y,z))),b),(c,d))`
Then the output tree will be: `(((a,(x,(y,z)<1>)<2>)<3>,b)<4>,(c,d)<5>)<6>`

#### `--multree` : Build MUL-trees

This option can be used with `-s`, `-h1`, and `-h2` to build MUL-trees from a standard species tree.

`python grampa.py -s species.tree -h1 "2" -h2 "4" -o multree.txt --multree`

If the input in `species.tree` is : (((a,(x,(y,z))),b),(c,d))
Then the output MUL-tree will be: ((((a,(x+,(y+,z+)<1>)<2>)<3>,b)<4>,(x*,(y*,z*)<5>)<6>)<7>,(c,d)<8>)<9>

#### `--tests` : Tests script

This option runs the tests.py script to ensure that all the other options in the program are working.

## License

This file is part of GRAMPA.

GRAMPA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

GRAMPA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 






