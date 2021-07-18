# GRAMPA
### Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis

## Authors
#### Gregg Thomas, S. Hussain Ather, Matthew Hahn

## About

### GRAMPA is a program to identify and place polyploidy events on a phylogeny and count duplications and losses in the presence of polyploidy.

##### What follows is a brief explanation of the options. Please see the project website for much more in depth info about the algorithm and examples of the program's usage.

#### https://gwct.github.io/grampa/

## Citation

#### Thomas GWC, Ather SH, Hahn MW. 2017. Gene-tree reconciliation with MUL-trees to resolve polyploid events. Systematic Biology. 66(6):1007-1018. https://doi.org/10.1093/sysbio/syx044

## Version History
#### This is version 1.3, released February 14, 2020
##### Note: as of this version I am recommending Python 3. Python 2 will no longer be supported going forward.

Change log:
* Added the `--numtrees` option to quickly calculate the number of MUL-trees possible for a given input.
* Added the `--version` option to quickly print out the version info.
* Implemented error checking for the psutil module when `--stats` is used.
* General upkeep.

###### Version 1.2 (April 16, 2017): Added the `-p` option so users can specify the number of processes GRAMPA should utilize.
###### Version 1.1 (March 22, 2017): Gene tree filtering step added and several other functions.
###### Version 1.0 (Summer, 2016): First release and implementation of MUL-tree reconciliation algorithm.


## Installation

Simply download the program and run it. You may want to add the GRAMPA folder to your $PATH variable for ease of use.
### The only dependency is Python 3 or higher

GRAMPA is also available through [bioconda](https://anaconda.org/bioconda/grampa). Thanks to [nathanweeks](https://github.com/nathanweeks) for setting that up!

## Usage

The first thing you should do when you try to run GRAMPA is make sure everything is working with some test files. You can do this easily by running the --tests command:

`python grampa.py --tests`

If all tests pass, then you're good to go! Basic usage in a real case would be:

`python grampa.py -s [species tree file] -g [gene trees file] -o [output directory]`

This would perform a full search for the optimal (lowest-scoring) MUL-tree on the input species tree. The other important options are `-h1` and `-h2`. Read more about them below and on the [project website](https://gwct.github.io/grampa/).


### Options

| Option | Description | 
| ------ | ----------- |
| -s | A file or string containing a bifurcating, rooted species tree in Newick format. This tree can either be singly-labeled or MUL. |
| -g | A file containing one or more bifurcating, rooted, Newick formatted gene trees. Gene trees with polytomies will be removed from the dataset. |
| -h1 | A space separated list of nodes to search as the polyploid clade. Only used with `-t s`. If nothing is entered all nodes will be considered. |
| -h2 | A space separated list of nodes to search as possible parental lineages for all nodes specified with `-h1`. If nothing is entered all possible nodes for the current `h1` will be considered. |
| -d | An option to specify whether to do reconciliations to MUL-trees only (0), the singly-labeled tree only (1), or both (2). Default: 2 |
| -c | The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 18 |
| -o | Output directory name. If the directory is not present, GRAMPA will created it for you. |
| -f | By default, all output files created by GRAMPA will have the prefix 'grampa_'. You can specify a different prefix with this option. |
| -v | Control the amount of output printed to the screen. Print all output (1) or just some log info (0). Default: 1 |
| -p | The number of processes GRAMPA should use |
| --multree | Set this flag if your input species tree is a MUL-tree. |
| --labeltree | The program will simply label your input species tree. |
| --numtrees | The program will simply count the number of possible MUL-trees given `-s`. `-h1`, and `-h2` may also be supplied.
| --buildmultrees | Build MUL-trees given `-s`, `-h1` and `-h2`. |
| --checknums | If this flag is entered, the program will just calculate the number of groups per gene tree and exit. No reconciliations will be done. |
| --maps | Output the node maps for each reconciliation in addition to the scores. The maps will be placed in the detailed output file. |
| --version | Print out version info and exit. |
| --tests | Run the tests script. |

