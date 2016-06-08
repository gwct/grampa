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

### Input

There are two main inputs for the program. 

1. A file containing a species tree (`-s`)
2. A file containing a list of gene trees (one tree per line). (`-g`)

The labels of the gene tree MUST be formatted such that they end with _[species label], where [species_label] corresponds to a tip label in the species tree.

### Output

GRAMPA creates two main output files: The main one you specify with -o gives the total reconciliation score for each MUL-tree considered, along with information about the trees. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.
The secondary output file (_det) is a detailed output describing the reconciliation scores from each gene tree to each MUL-tree.

### Options

| Option | Description | 
| ------ | ----------- |
| -s | A file containing a bifurcating species tree in newick format. This tree can either be standard or MUL |
| -t | This specifies the type of tree entered in `-s`. m: MUL-tree, s: standard tree. Default: s |
| -g | A file containing one ofr more newick formatted gene trees |
| -h1 | A comma separated list of nodes to search as the polyploid clade. Only used with `-t s` |
| -h2 | A comma separated list of nodes to search as possible parental lineages for all nodes specified with `-h1` |
| -c | The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 15 |
| -o | Output file name |
| --checknum | If this flag is entered, the program will just calculate the number of groups per gene tree. No reconciliations will be done |
| --lebeltree | The program will simply label your input species tree |
| --multree | Build MUL-trees given `-s`, `-h1`, and `-h2` |
| --tests | Run the tests script |
