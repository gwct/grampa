############################################################
# For GRAMPA site, 12.19
# This generates the file "readme.html"
############################################################

import sys, os
sys.path.append('..')
import lib.read_chunks as RC

######################
# HTML template
######################

html_template = """
<!doctype html>
    {head}

<body>
    {nav}

<div class="pure-g"><div class="pure-u-1" id="divider_row"></div></div>

	<div class="pure-g" id="main_row">
		<div class="pure-u-3-24" id="margin"></div>
		<div class="pure-u-18-24" id="main_col">
			<div id="main_content">
				<h1>README</h1>
				<h4>This section describes the program and its usage. For background about the algorithm see the <a href="readme.html">About</a> section.</h4>
				<div class="readme_divider"></div>
				<h2>GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis</h2>

				<div class="pure-g">
					<a name="install"></a>
					<div class="pure-u-1" id="jump_row">
						<div id="jump_container">
							<center>Jump to section:
							<a class="jump_link" href="#install">Installation and Usage</a>
							<a class="jump_link" href="#options">Options summary</a>
							<a class="jump_link" href="#detail_options">Options details</a>
							</center>
						</div>
					</div>
				</div>	

				<h2>Installation</h2>
				<div class="readme_divider"></div>
				<p>Clone or download the github repo: <a href="https://github.com/gwct/grampa" target="_blank">GRAMPA github</a></p>
				<p>The only dependency is Python 2.7 or higher. You may want to add the GRAMPA folder to your $PATH variable for ease of use!</p>

				<h2>Usage</h2>
				<div class="readme_divider"></div>
				<p>The first thing you should do when you try to run GRAMPA is make sure everything is working with some test files. You can do this easily 
					by running the <code class="cb">--tests</code> command:</p>

				<center><code class="cb">python grampa.py --tests</code></center>

				<p>If all tests pass, then you're good to go! Basic usage in a real case would be:</p>

				<center><code class="cb">python grampa.py -s [species tree file] -g [gene trees file] -o [output directory]</code></center>

				<p>This would perform a full search for the optimal (lowest scoring) MUL-tree on the input species tree.</p>

				<a name="input"></a><h3>Input</h3>
				<p>There are two main inputs for the program:</p>
				<ol>
					<li>A file or string containing a Newick formatted <b>rooted</b> species tree (<code class="cb">-s</code>). This can be a singly labeled tree 
						or a MUL-tree.</li>
					<li>A file containing one or more Newick formatted <b>rooted</b> gene trees (one tree per line) (<code class="cb">-g</code>).</li>
				</ol>

				<p>Important: the tip labels of the gene tree MUST be formatted such that they end with _[species label], where [species label] corresponds
					to a tip label in the species tree.</p>

				<a name="output"></a><h3>Output</h3>
				<p>All output files will be placed in the directory specified with <code class="cb">-o</code></p>
				<p>GRAMPA creates three main output files and a filtered tree file (if necessary).</p>
				<p>GRAMPA also creates a directory within the output directory called <b>groups_dir</b>. This just stores the gene tree groupings for each MUL-tree
					(in <a href="https://docs.python.org/2/library/pickle.html" target="_blank">pickled</a> format) so GRAMPA doesn't eat up a lot of RAM during
					reconciliations. This can be ignored/deleted</p>
				<ol>
					<li><h4>grampa_out.txt</h4>
						<p>This is the main output file and contains some log info for the current run and it gives the total reconciliation score for each MUL-tree
							considered. At the bottom of the file it will display the MUL-tree with the minimum reconciliation score.</p>

<pre><code># Tree #    H1 node H2 node Tree string Total score
ST          (((a,(x,(y,z)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,d)&lt;5&gt;)&lt;6&gt; 115
MT-1    &lt;2&gt; a   ((((a,(x*,(y*,z*)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,(x+,(y+,z+)&lt;4&gt;)&lt;5&gt;)&lt;6&gt;,b)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;  119
MT-2    &lt;2&gt; c   (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,((c,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,d)&lt;8&gt;)&lt;9&gt;  96
MT-3    &lt;2&gt; b   (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,(b,(x*,(y*,z*)&lt;4&gt;)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;  130
MT-4    &lt;2&gt; d   (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,(d,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;  76
MT-5    &lt;2&gt; &lt;3&gt; ((((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,(x*,(y*,z*)&lt;4&gt;)&lt;5&gt;)&lt;6&gt;,b)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;  119
MT-6    &lt;2&gt; &lt;2&gt; (((a,((x+,(y+,z+)&lt;1&gt;)&lt;2&gt;,(x*,(y*,z*)&lt;3&gt;)&lt;4&gt;)&lt;5&gt;)&lt;6&gt;,b)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;  145
MT-7    &lt;2&gt; &lt;6&gt; ((((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,d)&lt;5&gt;)&lt;6&gt;,(x*,(y*,z*)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;  132
MT-8    &lt;2&gt; &lt;4&gt; ((((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;  118
MT-9    &lt;2&gt; &lt;5&gt; (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,((c,d)&lt;5&gt;,(x*,(y*,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;  95
# ---------
The MUL-tree with the minimum parsimony score is MT-4:  (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,(d,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;
Score = 76
</code></pre>

						<ul>
							<li>The first line has the headers for the columns of the output table. Note that H1 node and H2 node are always nodes in the singly
								labeled tree.</li>
							<li>The total reconciliation score in the last column is the sum of all reconciliation scores for all gene trees for that MUL-tree.</li>
						</ul>
					</li>

					<li><h4>grampa_det.txt</h4>
						<p>The secondary output file contains detailed output describing the reconciliation scores from each gene tree to the lowest scoring MUL-tree.</p>

<pre><code># MT-4:(((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,(d,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;   H1 Node:&lt;2&gt; H2 Node:d
# GT/MT combo   # dups  # losses    Total score
GT-1 to MT-4    1   1   2
GT-2 to MT-4    1   2   3
GT-3 to MT-4    1   3   4
GT-4 to MT-4    1   2   3
.
.
.
GT-25 to MT-4   0   1   1
# Gene trees with multiple maps:    1
# Total parsimony score for MT-4: 76
</code></pre>
					
						<ul>
							<li>The first line contains info about the MUL-tree and the second line contains the headers for the rest of the table</li>
							<li>Note that the lowest score for some GT/MT combos can have multiple maps. In these cases, we report all possible scores</li>
						</ul>
					</li>

					<li><h4>grampa_checknums.txt</h4>
						<p>GRAMPA must calculate how many combinations of maps there are for each gene-tree/MUL-tree pair and filter out those that are over
							the group cap in any combo before any reconciliations can be done. This filtering ensures that all MUL-trees are reconciled to
							the same set of gene-trees. The number of groups for each gene-tree/MUL-tree is recorded in this file.</p>

<pre><code># GT/MT combo   # Groups    # Fixed # Combinations
# MT-1:((((a,(x*,(y*,z*)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,(x+,(y+,z+)&lt;4&gt;)&lt;5&gt;)&lt;6&gt;,b)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;   H1 Node:&lt;2&gt; H2 Node:a
GT-1 to MT-1    2   0   4
GT-2 to MT-1    1   1   2
GT-3 to MT-1    1   1   2
GT-4 to MT-1    1   1   2
.
.
.
GT-25 to MT-9   0   2   1
</code></pre>

						<ul>
							<li>The first line in the file contains the table headers. For each MUL-tree there is also a line giving som info about it.</li>
							<li>If a gene tree is over the cap specified with <code class="cb">-c X</code> then the message <code class="cb">Number of 
								groups over group cap (-c set to X) -- Filtering.</code> will also be displayed on the line.</li>
						</ul>
					</li>
				</ol>

				<div class="pure-g">
					<a name="options"></a>
					<div class="pure-u-1" id="jump_row">
						<div id="jump_container">
							<center>Jump to section:
							<a class="jump_link" href="#install">Installation and Usage</a>
							<a class="jump_link" href="#options">Options summary</a>
							<a class="jump_link" href="#detail_options">Options details</a>
							</center>
						</div>
					</div>
				</div>	

				<h2>Options Table</h2>
				<div class="readme_divider"></div></br>
					<table class="pure-table pure-table-bordered pure-table-striped">
						<thead><tr><th>Option</th><th>Description</th></tr></thead>
						<tbody>
							<tr>
								<td>-s</td>
								<td>A file or string containing a bifurcating, rooted species tree in Newick format. This tree can either be singly-labeled or MUL.</td>
							</tr>
							<tr>
								<td>-g</td>
								<td>A file containing one or more rooted, Newick formatted gene trees.</td>
							</tr>
							<tr>
								<td>-h1</td>
								<td>A space separated list of nodes to search as the polyploid clade. If nothing is entered all nodes will be considered.</td>
							</tr>
							<tr>
								<td>-h2</td>
								<td>A space separated list of nodes to search as possible parental lineages for all nodes specified with <code class="cb">-h1</code>. 
									If nothing is entered all possible nodes for the current <code class="cb">h1</code> will be considered.</td>
							</tr>
							<tr>
								<td>-d</td>
								<td>An option to specify whether to do reconciliations to MUL-trees only (0), the singly-labeled tree only (1), or both (2). Default: 2</td>
							</tr>
							<tr>
								<td>-c</td>
								<td>The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 18</td>
							</tr>
							<tr>
								<td>-o</td>
								<td>Output directory name. If the directory is not present, GRAMPA will created it for you.</td>
							</tr>
							<tr>
								<td>-f</td>
								<td>By default, all output files created by GRAMPA will have the prefix 'grampa_'. You can specify a different prefix with this option.</td>
							</tr>
							<tr>
								<td>-v</td>
								<td>Control the amount of output printed to the screen. Print all output (1) or just some log info (0). Default: 1</td>
							</tr>
							<tr>
								<td>-p</td>
								<td>The number of processes GRAMPA should use for reconciliations.</td>
							</tr>
							<tr>
								<td>--multree</td>
								<td>Set this flag if your input species tree is a MUL-tree.</td>
							</tr>
							<tr>
								<td>--labeltree</td>
								<td>The program will simply label your input species tree.</td>
							</tr>
							<tr>
								<td>--buildmultrees</td>
								<td>Build MUL-trees given <code class="cb">-s</code>, <code class="cb">-h1</code>, and <code class="cb">-h2</code>.</td>
							</tr>
							<tr>
								<td>--checknums</td>
								<td>If this flag is entered, the program will just calculate the number of groups per gene tree and exit. No reconciliations will be done.</td>
							</tr>
							<tr>
								<td>--maps</td>
								<td>Output the node maps for each reconciliation in addition to the scores. The maps will be placed in the detailed output file.</td>
							</tr>
							<tr>
								<td>--tests</td>
								<td>Run the tests script</td>
							</tr>
						</tbody>
					</table>
					</br>
				<div class="pure-g">
					<a name="detail_options"></a>
					<div class="pure-u-1" id="jump_row">
						<div id="jump_container">
							<center>Jump to section:
							<a class="jump_link" href="#install">Installation and Usage</a>
							<a class="jump_link" href="#options">Options summary</a>
							<a class="jump_link" href="#detail_optoins">Options details</a>
							</center>
						</div>
					</div>
				</div>

				<h2>Detailed options</h2>
				<div class="readme_divider"></div></br>
					<a name="-s"></a><h3><code class="cb">-s</code> : A rooted, Newick formatted species tree. This tree can be singly-labeled or MUL.</h3>
						<ul>
							<p>The tree can be in a file, in which case you enter the file name here, or you can simply paste the tree string into the command line.</p>
							<p>Entering a singly-labeled tree means you wish to search for the most parsimonious polyploidy scenario. GRAMPA will build MUL-trees based
								on this singly-labeled tree and calculate reconciliation scores. You can specify the range of MUL-trees to build with the <code class="cb">
								-h1</code> and <code class="cb">-h2</code> options.</p>
							<p>Example singly-labeled species tree:</p>
							<center><code class="cb">(((a,(x,(y,z))),b),(c,d))</code></center>
							<p>Entering a MUL-tree is the equivalent of entering a singly-labeled tree and specifying a single H1 and single H2 node. It represents a 
								single scenario of polyploidy and can be used if you wish to count the number of duplications and losses on gene trees given that
								scenario.</p>
							<p><b>NOTE: If a MUL-tree is entered, the <code class="cb">--multree</code> flag must be set.</b></p>
							<p>Example MUL-tree:</p>
							<center><code class="cb">((((a,(x,(y,z))),b),(x,(y,z))),(c,d))</code></center>
						</ul>

					<a name="-g"></a><h3><code class="cb">-g</code> : A file containing newick formatted gene trees.</h3>
						<ul>
							<p>This file should contain one or more Newick formatted gene trees, with one tree per line in the file.</p>
							<p><b>The tip labels in the gene trees must end with _[species label]</b> where [species label] matches a tip label in the species tree
								This is necessary so GRAMPA can initialize the mappings correctly.</p>
							<p>Alternatively, if you wish to reconcile to only a single gene tree, you can simply paste the tree string into the command line.</p>
						</ul>

					<a name="h1h2"></a><h3><code class="cb">-h1</code> and <code class="cb">-h2</code> : GRAMPA's search parameters.</h3>
						<ul>
							<p>H1 and H2 are nodes in the singly-labeled species tree that define how to build a MUL-tree. H1 is the node that represents the
								polyploid clade. The subtree rooted at H1 and the branch that H1 subtends will be copied onto the branch that H2 subtends:</p>
							<img class="pure-img" id="logo_main" src="img/mul_tree.png">
							<p>In the above example, H1 is node 2 and H2 is node 5 in the singly-labeled tree. This leads to the MUL-tree on the right.</p>
							<p>H1 and H2 can be input in 2 different, equivalent ways:</p>
							<center><code class="cb">-h1 "2" -h2 "5"</code> and <code class="cb">-h1 "x,y,z" -h2 "c,d"</code></center>
							<p>The first way relies on internal node labels. To label your species tree, use the <code class="cb">--labeltree</code> option.
								</br><b>IMPORTANT: For now, only use node labels as specified by <code class="cb">--labeltree</code>. Custom labels will not work.</b></p>
							<p>The second way uses a list of the species that define that node. <b>Species must be comma delimited.</b></p>
							<p>H2 cannot be located below H1 in the species tree! If this occurs, GRAMPA will just tell you that it's not possible and move on.</p>
							<p>Multiple H1 and H2 nodes can be entered as a space delimited list:</p>
							<center><code class="cb">-h1 "2 3" -h2 "5 6"</code> and <code class="cb">-h1 "x,y,z a,x,y,z" -h2 "c,d a,b,c,d,x,y,z"</code> are equivalent.</center>
							<p>Entering this means that GRAMPA will first set H1 as node 2 and try both nodes 5 and 6 as H2. Then H1 will be set to node 3 and will try nodes 5 
								and 6 as H2.</p>
							<h4>If <code class="cb">-h1</code> and <code class="cb">-h2</code> are not specified, GRAMPA will try all possible node combinations
								of H1 and H2!</h4>
						</ul>

					<a name="-d"></a><h3><code class="cb">-d</code> : Reconciliation type option</h3>
						<ul>
							<p>GRAMPA can do reconciliations to singly-labeled and MUL-trees. If you know a polyploidy event has taken place, you may wish to only 
								reconcile to MUL-trees. However, if you are trying to identify a new polyploidy event, the scores of all MUL-trees considered must 
								be compared to the score of the singly-labeled tree, which represents a scenario of no polyploidy.</p>
							<center>
							<table class="pure-table pure-table-bordered pure-table-striped">
								<thead><tr><th>Input</th><th>Setting</th></tr></thead>
								<tbody>
									<tr><td><code class="cb">-d 0</code></td><td>Reconcile to MUL-trees only</td></tr>
									<tr><td><code class="cb">-d 1</code></td><td>Reconcile to the singly-labeled tree only</td></tr>
									<tr><td><code class="cb">-d 2</code></td><td>Reconcile to both the singly-labeled tree and MUL-trees</td></tr>
								</tbody>
							</table>
							</center>
							<p>Setting <code class="cb">-d 1</code> also means you can use GRAMPA to count duplications and losses in the absence of polyploidy, 
								like any other reconciliation program!</p>
						</ul>

					<a name="-c"></a><h3><code class="cb">-c</code> : The group cap</h3>
						<ul>
							<p>GRAMPA uses the standard LCA reconciliation algorithm on MUL-trees, meaning that some genes have more than one possible mapping.
								We get around this by trying ALL possible initial mappings and picking the one with the lowest score. This works, but also means
								our program has an exponential runtime based on the number of genes from polyploid species in any given gene tree. We get around
								this in several ways by collapsing and fixing groups (see paper), but there can still be lots of groups. This parameter sets the
								maximum number of groups to consider for any gene tree. If a gene tree has more than this number of groups, it will be skipped.</p>
							<p>Default is 8 groups, with a max setting of 18.</p>
						</ul>

					<a name="-o"></a><h3><code class="cb">-o</code> : Output directory</h3>
						<ul>
							<p>Grampa creates several output files, so it is easiest just to place them all in a single directory. That directory can be specified
								with this option, and will be created for you if it doesn't exist. If this option is not specified, the default output directory
								is "grampa_[date]-[time]".</p>
						</ul>

					<a name="-f"></a><h3><code class="cb">-f</code> : Output file prefix</h3>
						<ul>
							<p>By default, all output files created by GRAMPA will have the prefix 'grampa_'. You can specify a different prefix with this option.
								For example, a run with <code class="cb">-p test</code> will generate the following output files, all within the output directory:</p>

							<center><code class="cb">test_out.txt, test_det.txt, test_checknums.txt</code></center>
						</ul>

					<a name="--multree"></a><h3><code class="cb">--multree</code> : Input MUL-tree flag</h3>
						<ul>
							<p>GRAMPA can accept both singly-labeled and MUL-trees as input. If your input species tree (<code class="cb">-s</code>) is a MUL-tree,
								you must set this flag so GRAMPA knows to read it as a MUL-tree. A MUL-tree represents a single possible polyploid scenario and it is
								equivalent to entering a singly-labeled tree with a single H1 and H2 node specified.</p>
						</ul>

					<a name="--labeltree"></a><h3><code class="cb">--labeltree</code> : Species tree labeling</h3>
						<ul>
							<p>This option can be used in conjunction with <code class="cb">-s</code> to simply add internal node labels to a species tree and print
								it to the screen. For example, if the file <code class="cb">species.tree</code> contains the following tree:</p>

							<center><code class="cb">(((a,(x,(y,z))),b),(c,d))</code></center>

							<p>Then the command:</p>

							<center><code class="cb">python grampa.py -s species.tree --labeltree</code></center>

							<p>Will simply print this to the screen as output:</p>

							<center><code class="cb">(((a,(x,(y,z)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,d)&lt;5&gt;)&lt;6&gt;</code></center>
						</ul>

					<a name="--buildmultrees"></a><h3><code class="cb">--buildmultrees</code> : Building MUL-trees</h3>
						<ul>
							<p>This option can be used with <code class="cb">-s</code>, <code class="cb">-h1</code>, and <code class="cb">-h2</code> to build
								MUL-trees from a standard species tree. For example, if the file <code class="cb">species.tree</code> contains the following tree:</p>

							<center><code class="cb">(((a,(x,(y,z))),b),(c,d))</code></center>

							<p>Then the command:</p>

							<center><code class="cb">python grampa.py -s species.tree -h1 "2" -h2 "4" -o multree_ex --buildmultrees</code></center>

							<p>Will yield the following output in the main output file (multree_ex/grampa_out.txt):</p>

							<center><code class="cb">((((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;</code></center>
						</ul>

					<a name="--checknums"></a><h3><code class="cb">--checknums</code> : Group counting</h3>
						<ul>
							<p>With this set, the program will run normally with the specified options, except no reconciliations will be done. Instead, only the
								checknums output file will be created and will contain the number of polyploid groups for each gene tree. Use this to decide the
								best setting for <code class="cb">-c</code>.</p>
						</ul>

					<a name="--maps"></a><h3><code class="cb">--maps</code> : Output node mappings</h3>
						<ul>
							<p>Set this option to output the LCA node mappings along with the reconciliation scores to the detailed output file. This adds a column
							to the _det.txt output file with that gene tree with the nodes re-labeled to include the maps, dups, and losses along that branch.</p>

<pre><code># GT/MT combo	# dups	# losses	Total score	Maps
GT-1 to MT-4	1	1	2	(((1_a[a+-0],((1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;3&gt;[&lt;1&gt;+-0])&lt;4&gt;[&lt;2&gt;+-0])&lt;5&gt;[&lt;2&gt;+-1])&lt;6&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;7&gt;[&lt;4&gt;+-0],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
GT-2 to MT-4	1	2	3	((((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;4&gt;[&lt;1&gt;+-0])&lt;5&gt;[&lt;2&gt;+-0])&lt;6&gt;[&lt;3&gt;+-1],1_b[b+-0])&lt;7&gt;[&lt;4&gt;+-0],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
GT-3 to MT-4	1	3	4	((((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;4&gt;[&lt;4&gt;+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;5&gt;[&lt;1&gt;+-0])&lt;6&gt;[&lt;2&gt;+-0])&lt;7&gt;[&lt;4&gt;+-1],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
GT-4 to MT-4	1	2	3	((((1_a[a+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;4&gt;[&lt;1&gt;+-0])&lt;5&gt;[&lt;2&gt;+-0])&lt;6&gt;[&lt;3&gt;+-1],1_b[b+-0])&lt;7&gt;[&lt;4&gt;+-0],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]

.
.
.
GT-25 to MT-4	0	1	1	(((1_a[a+-0],1_x[x+-0])&lt;1&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;2&gt;[&lt;4&gt;+-0],(1_c[c+-0],(1_d[d+-0],(2_x[x*-0],(1_y[y*-0],1_z[z*-0])&lt;3&gt;[&lt;5&gt;+-0])&lt;4&gt;[&lt;6&gt;+-0])&lt;5&gt;[&lt;7&gt;+-0])&lt;6&gt;[&lt;8&gt;+-0])&lt;7&gt;[&lt;9&gt;+-0]
# Gene trees with multiple maps:    1
# Total parsimony score for MT-4: 76
</code></pre>
							
							<p>The last column is the gene tree with the nodes relabeled as:</p>
							<center><code class="cb">Node[Map-Dups]</code></center>
							<p>where Dups is 1 if the node is a duplication node and 0 if not.</p>
							<p>These trees can be rendered with a tree viewer such as <a href="http://doua.prabi.fr/software/seaview" target="_blank">SeaView</a>
								or <a href="http://tree.bio.ed.ac.uk/software/figtree/" target="_blank">FigTree</a>.</p>
						</ul>
			</div>
		</div>
		<div class="pure-u-3-24" id="margin"></div>
	</div>

    {footer}
</body>
"""

######################
# Main block
######################
pagefile = "readme.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - README"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));