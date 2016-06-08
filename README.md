# GRAMPA
### Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis

## Author
#### Gregg Thomas

## About

#### GRAMPA is a program to date polyploidy events and count duplications and losses in the presence of polyploidy.

#### This is version Beta 1.0
#### The only dependency is Python 2.7 or higher

## Installation

Simply download the program and run it. You may want to add the GRAMPA folder to your $PATH variable for ease of use.

## Usage

The first thing you should do when you first try to run GRAMPA is make sure everything is working with some test files. You can do this easily by running the --tests command:

`python grampa.py --tests`

If all tests pass, then you're good to go! Basic usage in a real case would be:

`python grampa.py -s [species tree file] -g [gene trees file] -o [output file name]`

This would perform a full search for the optimal MUL-tree on the species tree.

### Input

There are two main inputs for the program. 

1. A file containing a species tree (`-s`)
2. A file containing a list of gene trees (one tree per line). (`-g`)

The labels of the gene tree MUST be formatted such that they end with _[species label], where [species_label] corresponds to a tip label in the species tree.

### Output

GRAMPA creates two main output files, both specified with `-o`: 

1. The main output gives the total reconciliation score for each MUL-tree considered, along with information about the trees. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.

2. The secondary output file is a detailed output describing the reconciliation scores from each gene tree to each MUL-tree. It has the same name as the file specified with `-o`, but with the suffix '_det' added to the end.

### Options

| Option | Description | 
| ------ | ----------- |
| -s | A file containing a bifurcating species tree in newick format. This tree can either be standard or MUL |
| -t | This specifies the type of tree entered in `-s`. m: MUL-tree, s: standard tree. Default: s |
| -g | A file containing one or more newick formatted gene trees |
| -h1 | A space separated list of nodes to search as the polyploid clade. Only used with `-t s`. If nothing is entered all nodes will be considered |
| -h2 | A space separated list of nodes to search as possible parental lineages for all nodes specified with `-h1`. If nothing is entered all possible nodes for the current `h1` will be considered |
| -c | The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 15 |
| -o | Output file name |
| --checknum | If this flag is entered, the program will just calculate the number of groups per gene tree. No reconciliations will be done |
| --labeltree | The program will simply label your input species tree |
| --multree | Build MUL-trees given `-s`, `-h1`, and `-h2` |
| --tests | Run the tests script |

#### Detailed options

###### `-s` : A file containing a newick formatted species tree. This tree can be either standard or MUL. 

Entering a standard tree means you wish to search for the most parsimonious polyploidy scenario. GRAMPA will build MUL-trees based on this standard tree and calculate reconciliation scores. You can specify the range of MUL-trees to build with the `-h1` and `-h2` options.

Example standard tree: `(((a,(x,(y,z))),b),(c,d))`

Entering a MUL-tree is the equivalent of entering a standard tree and specifying a single H1 and single H2 node. It represents a single scenario of polyploidy and should be used if you wish to count the number of duplications and losses on gene trees given that scenario. **NOTE: If a MUL-tree is entered, you MUST set `-t m`.**

Example MUL-tree: `((((a,(x,(y,z))),b),(x,(y,z))),(c,d))`

###### `-t` : Input species tree type.

`-t s` : The tree specified with `-s` is standard. This is the Default setting.
`-t m` : The tree specified with `-s` is a MUL-tree.

###### `-g` : A file containing newick formatted gene trees.

This file should contain one or more gene trees, with one tree per line in the file. 

**The labels of the gene trees MUST end with `_[species label]`**

Where `[species labe]` matches a tip label in the species tree. This is necessary so that we can initialize the mappings correctly.

###### `-h1` and `-h2` : GRAMPA's search parameters

H1 and H2 are nodes in the species tree that define how to build a MUL-tree:

-> ![test](https://github.com/gwct/grampa/blob/master/test2.png) <-





