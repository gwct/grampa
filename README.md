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

<p align="center"> `python grampa.py --tests` </p>

### Input

The main inputs of the program are a file containing a species tree and a file containing a list of gene trees (one tree per line).
The labels of the gene tree MUST be formatted such that they end with _[species label], where [species_label] corresponds to a tip label in the species tree.

### Output

GRAMPA creates two main output files: The main one you specify with -o gives the total reconciliation score for each MUL-tree considered, along with information about the trees. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.
The secondary output file (_det) is a detailed output describing the reconciliation scores from each gene tree to each MUL-tree.
