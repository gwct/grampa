# GRAMPA
### Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis

## Author
#### Gregg Thomas

## About

#### This is version Beta 1.0
#### The only dependency is Python 2.7

## Usage

### Input

The main inputs of the program are a file containing a species tree and a file containing a list of gene trees (one tree per line).
The labels of the gene tree MUST be formatted such that they end with _[species label], where [species_label] corresponds to a tip label in the species tree.

### Output

GRAMPA creates two main output files: The main one you specify with -o gives the total reconciliation score for each MUL-tree considered, along with information about the trees. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.
The secondary output file (_det) is a detailed output describing the reconciliation scores from each gene tree to each MUL-tree.
