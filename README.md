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
#### This is version 1.4.0, released March 24, 2023

Change log:
* Re-factor of much of the code that handles inputs and outputs.
* Reformatted output files to be easier to handle when post-processing (see the [website README](https://gwct.github.io/grampa/readme.html) for more info).
* Removed `-d` in favor of `--st-only` and `--no-st`

###### Version 1.3 (February 14, 2020): Added the `--numtrees` and `--version` options. Dropped support for Python 2.
###### Version 1.2 (April 16, 2017): Added the `-p` option so users can specify the number of processes GRAMPA should utilize.
###### Version 1.1 (March 22, 2017): Gene tree filtering step added and several other functions.
###### Version 1.0 (Summer, 2016): First release and implementation of MUL-tree reconciliation algorithm.


## Installation

### Install from [bioconda](https://anaconda.org/bioconda/grampa)

```
mamba install grampa
```

*Note that if you prefer `conda`, simply replace `mamba` with `conda` in the command above.*

Once installed, the program can be called as `grampa` or `grampa.py`

*Thanks to [nathanweeks](https://github.com/nathanweeks) for setting up the bioconda recipe!*

### Install from source

Simply download the [latest release](https://github.com/gwct/grampa/releases/latest), decompress, and run it. GRAMPA is a Python script with no external dependencies so it should work out of the box.

*Note that if installed from source, you will have to invoke Python explicitly and provide the path to the GRAMPA script when you run it: `python /path/to/grampa.py`. You may want to add the GRAMPA folder to your $PATH variable for ease of use.*

## Usage

*The following example commands assume an install from bioconda and use `grampa` to call the program. If you installed from source, see above for information on how to call the program.*

If necessary, you can run GRAMPA with some test files. You can do this easily by running the --tests command:

`grampa --tests`

If all tests pass, then you're good to go! Basic usage in a real case would be:

`grampa -s [species tree file] -g [gene trees file] -o [output directory]`

This would perform a full search for the optimal (lowest-scoring) MUL-tree on the input species tree. The other important options are `-h1` and `-h2`. Read more about them below and on the [project website](https://gwct.github.io/grampa/).


### Options

| Option | Description | 
| ------ | ----------- |
| -s | A file or string containing a bifurcating, rooted species tree in Newick format. This tree can either be singly-labeled or MUL. |
| -g | A file containing one or more bifurcating, rooted, Newick formatted gene trees. Gene trees with polytomies will be removed from the dataset. |
| -h1 | A space separated list of nodes to search as the polyploid clade. Cannot be used with used with `--multree`. If nothing is entered all nodes will be considered. |
| -h2 | A space separated list of nodes to search as possible parental lineages for all nodes specified with `-h1`.  Cannot be used with used with `--multree`. If nothing is entered all possible nodes for the current `h1` will be considered. |
| -c | The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 18 |
| -o | Output directory name. If the directory is not present, GRAMPA will created it for you. |
| -f | By default, all output files created by GRAMPA will have the prefix 'grampa'. You can specify a different prefix with this option. |
| -v | Control the amount of output printed to the screen. 0: print nothing. 1: print only some info at the start. 2: print all log info to screen. 3 (default): print all log info to the screen as well as progress updates for certain steps. |
| -p | The number of processes GRAMPA should use. |
| --multree | Set this flag if your input species tree is a MUL-tree. |
| --labeltree | The program will simply label your input species tree. |
| --numtrees | The program will simply count the number of possible MUL-trees given `-s`. `-h1`, and `-h2` may also be supplied.
| --buildmultrees | Build MUL-trees given `-s`, `-h1` and `-h2`. |
| --checknums | If this flag is entered, the program will just calculate the number of groups per gene tree and exit. No reconciliations will be done. |
| --st-only | Only perform reconciliations to the singly-labeled tree input with `-s`. |
| --no-st | Skip performing reconciliations to the singly-labeled tree input with `-s` and only do reconciliations to MUL-trees. |
| --maps | Output the node maps for each reconciliation in addition to the scores. The maps will be placed in the detailed output file. |
| --version | Print out version info and exit. |
| --tests | Run the tests script. |

