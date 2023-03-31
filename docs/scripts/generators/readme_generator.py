############################################################
# For GRAMPA site, 12.19
# This generates the file "readme.html"
############################################################

import sys, os
sys.path.append(os.path.abspath('../lib/'))
import read_chunks as RC

######################
# HTML template
######################

html_template = """
<!doctype html>
    {head}

<body>
    {nav}

	<div class="row" id="main_row">
		<div class="col-3-24" id="margin"></div>
		<div class="col-18-24" id="main_col">
			<div id="main_content">
				<h1>README</h1>

				<h4>This section describes the program and its usage. For background about the algorithm see the <a href="readme.html">About</a> section.</h4>

				<div class="readme_divider"></div>
				<h2>GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis</h2>

				<div class="row">
					<a name="install"></a>
					<div class="col-1" id="jump_row">
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

				<p>The only dependency is Python 3 or higher. You may want to add the GRAMPA folder to your $PATH variable for ease of use!</p>

				<h2>Usage</h2>

				<div class="readme_divider"></div>

				<p>The first thing you should do when you try to run GRAMPA is make sure everything is working with some test files. You can do this easily 
					by running the <code class="inline">--tests</code> command:</p>

				<center><code class="cb">python grampa.py --tests</code></center>

				<p>If all tests pass, then you're good to go! Basic usage in a real case would be:</p>

				<center><code class="cb">python grampa.py -s [species tree file] -g [gene trees file] -o [output directory]</code></center>

				<p>This would perform a full search for the optimal (lowest scoring) MUL-tree on the input species tree.</p>

				<a name="input"></a><h3>Input</h3>

				<p>There are two main inputs for the program:</p>

				<ol>
					<li>A file or string containing a Newick formatted <b>rooted</b> species tree (<code class="inline">-s</code>). This can be a singly labeled tree 
						or a MUL-tree.</li>
					<li>A file containing one or more Newick formatted <b>rooted</b> gene trees (one tree per line) (<code class="inline">-g</code>).</li>
				</ol>

				<p>Important: the tip labels of the gene tree MUST be formatted such that they end with _[species label], where [species label] corresponds
					to a tip label in the species tree.</p>

				<a name="output"></a><h3>Output</h3>

				<p>All output files will be placed in the directory specified with <code class="inline">-o</code></p>

				<p>GRAMPA creates four output files, a log file, and a filtered tree file (if necessary).</p>

				<p>GRAMPA also creates a directory within the output directory called <b>groups_dir</b>. This just stores the gene tree groupings for each MUL-tree
					(in <a href="https://docs.python.org/2/library/pickle.html" target="_blank">pickled</a> format) so GRAMPA doesn't eat up a lot of RAM during
					reconciliations. This can be ignored/deleted</p>

				<ol>
					<li><h4>grampa-scores.txt</h4>
						<p>This is the main output file and contains the total reconciliation score for each MUL-tree considered, sorted in ascending order. </p>

<!--
<pre><code>mul.tree        h1.node h2.node score   labeled.tree
106     &lt;2&gt;     d       76      (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,(d,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;
114     &lt;3&gt;     d       86      (((a+,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,(d,(a*,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;)&lt;10&gt;
110     &lt;2&gt;     &lt;5&gt;     95      (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,((c,d)&lt;5&gt;,(x*,(y*,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;
120     &lt;4&gt;     d       95      (((a+,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b+)&lt;4&gt;,(c,(d,((a*,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,b*)&lt;8&gt;)&lt;9&gt;)&lt;10&gt;)&lt;11&gt;
105     &lt;2&gt;     c       96      (((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,((c,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,d)&lt;8&gt;)&lt;9&gt;
117     &lt;3&gt;     &lt;5&gt;     104     (((a+,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,((c,d)&lt;5&gt;,(a*,(x*,(y*,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;)&lt;10&gt;
113     &lt;3&gt;     c       106     (((a+,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,((c,(a*,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,d)&lt;9&gt;)&lt;10&gt;
122     &lt;4&gt;     &lt;5&gt;     113     (((a+,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b+)&lt;4&gt;,((c,d)&lt;5&gt;,((a*,(x*,(y*,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,b*)&lt;9&gt;)&lt;10&gt;)&lt;11&gt;
123     &lt;4&gt;     &lt;6&gt;     113     ((((a+,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b+)&lt;4&gt;,(c,d)&lt;5&gt;)&lt;6&gt;,((a*,(x*,(y*,z*)&lt;7&gt;)&lt;8&gt;)&lt;9&gt;,b*)&lt;10&gt;)&lt;11&gt;
0       NA      NA      115     (((a,(x,(y,z)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,d)&lt;5&gt;)&lt;6&gt;
</code></pre>
-->

						<p>
							The first line of this file contains the headers, defined as follows for each subsequent row:
                        </p>
                        
						<div class="table-container">
							<table class="grid-table">
								<thead>
                                	<tr>
										<th>mul.tree</th>
										<th>h1.node</th>
										<th>h2.node</th>
										<th>score</th>
										<th>mul.tree</th>
                                	</tr>
                                </thead>
								<tbody>
									<tr>
										<td>The ID of the MUL-tree</td>
                                        <td>The H1 node in the species tree for the current MUL-tree</td>
                                        <td>The H2 node in the species tree for the current MUL-tree</td>
                                        <td>The total parsimony score for the current MUL-tree</td>
                                        <td>The Newick formatted tree string for the MUL-tree, with hybrid clades indicated with *</td>
                                    </tr>
								</tbody>		
							</table>
						</div>

						<h3>
							Please note that the input singly-labeled species tree always has the ID of 0
                        </h3>

					<li><h4>grampa-detailed.txt</h4>
						<p>The secondary output file contains detailed output describing the reconciliation scores from each gene tree to the lowest scoring MUL-tree.</p>

<!--                        
<pre><code>mul.tree	gene.tree	dups	losses	total.score
106	1	1	1	2
106	2	1	2	3
106	3	1	3	4
106	4	1	2	3
106	5	1	3	4
106	6	1	3	4
106	7	0	0	0
106	8	1	3	4
</code></pre>
-->

						<p>
							The first line of this file contains the headers, defined as follows for each subsequent row:
                        </p>

						<div class="table-container">
							<table class="grid-table">
								<thead>
                                	<tr>
										<th>mul.tree</th>
										<th>gene.tree</th>
										<th>dups</th>
										<th>losses</th>
										<th>total score</th>
                                	</tr>
                                </thead>
								<tbody>
									<tr>
										<td>The ID of the MUL-tree</td>
                                        <td>The ID of the gene tree being reconciled to the MUL-tree</td>
                                        <td>The number of duplications on this gene tree given this MUL-tree</td>
                                        <td>The number of losses on this gene tree given this MUL-tree</td>
                                        <td>The sum of dups and losses for this gene tree and MUL-tree</td>
                                    </tr>
								</tbody>	
							</table>
						</div>

						<h3>
                        	Note that the lowest score for some GT/MT combos can have multiple maps. In these cases, we report all possible scores.
                        </h3>

                    <li><h4>grampa-dup-counts.txt</h4>
						<p>
							For the <b>6 lowest scoring MUL-trees</b>, GRAMPA counts the number of duplications along each branch in the MUL-tree summed
                            over all gene trees.
                        </p>
                        
						<p>
							The first line of this file contains the headers, defined as follows for each subsequent row:
                        </p>

						<div class="table-container">
							<table class="grid-table">
								<thead>
                                	<tr>
										<th>mul.tree</th>
										<th>node</th>
										<th>dups</th>
                                	</tr>
                                </thead>
								<tbody>
									<tr>
										<td>The ID of the MUL-tree</td>
                                        <td>The node ID for the current MUL-tree</td>
                                        <td>The total number of duplications over all gene tres along the branch above the node in the MUL-tree</td>
                                    </tr>
								</tbody>	
							</table>
						</div>                       
                        
					<li><h4>grampa-checknums.txt</h4>
						<p>GRAMPA must calculate how many combinations of maps there are for each gene-tree/MUL-tree pair and filter out those that are over
							the group cap in any combo before any reconciliations can be done. This filtering ensures that all MUL-trees are reconciled to
							the same set of gene-trees. The number of groups for each gene-tree/MUL-tree is recorded in this file.</p>

<!--
<pre><code>mul.tree        gene.tree       groups  fixed   combinations    over.cap.filtered
# MT-1:((((a+,a*)<1>,(x,(y,z)<2>)<3>)<4>,b)<5>,(c,d)<6>)<7>     H1 Node:a       H2 Node:a
1       1       1       0       2       N
1       2       1       0       2       N
1       3       1       0       2       N
1       4       1       0       2       N
1       5       1       0       2       N
1       6       1       0       2       N
</code></pre>
-->

						<p>
							The first line of this file contains the headers, defined as follows for each subsequent row:
                        </p>

						<div class="table-container">
							<table class="grid-table">
								<thead>
                                	<tr>
										<th>mul.tree</th>
										<th>gene.tree</th>
										<th>groups</th>
										<th>fixed</th>
										<th>combinations</th>
                                        <th>over.cap.filtered</th>
                                	</tr>
                                </thead>
								<tbody>
									<tr>
										<td>The ID of the MUL-tree</td>
                                        <td>The ID of the gene tree to be reconciled to the MUL-tree</td>
                                        <td>The number of distinct hybrid clades in the gene tree</td>
                                        <td>The number of hybrid clades in the gene tree that also group with a sister species from the singly-labeled tree</td>
                                        <td>The total number of mappings to try for the gene tree with this MUL-tree</td>
                                        <td>Either Y or N to indicate whether the number of groups exceeds the number set with -c</td>
                                    </tr>
								</tbody>	
							</table>
						</div>   
                        
                    <li><h4>grampa-trees-filtered.txt</h4>
						<p>
							A text file with the gene trees used for this GRAMPA run, after filtering by the group cap. One tree per line.
                        </p>
                    
                    <li><h4>grampa.log</h4>
						<p>
							A log file containing run time information and a summary of the lowest scoring MUL-tree.
                        </p>
				</ol>

				<div class="row">
					<a name="options"></a>
					<div class="col-1" id="jump_row">
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


					<div class="table-container">
						<table class="stripe-table">
							<thead>
								<tr>
									<th>Option</th>
									<th>Description</th>
								</tr>
							</thead>

							<tbody>
								<tr>
									<td><code class="inline">-s</code></td>
									<td>A file or string containing a bifurcating, rooted species tree in Newick format. This tree can either be singly-labeled or MUL.</td>
								</tr>
								<tr>
									<td><code class="inline">-g</code></td>
									<td>A file containing one or more bifurcating, rooted, Newick formatted gene trees. Gene trees with polytomies are currently not supported
										and will be automatically filtered from the analysis.</td>
								</tr>
								<tr>
									<td><code class="inline">-h1</code></td>
									<td>A space separated list of nodes to search as the polyploid clade. If nothing is entered all nodes will be considered.</td>
								</tr>
								<tr>
									<td><code class="inline">-h2</code></td>
									<td>A space separated list of nodes to search as possible parental lineages for all nodes specified with <code class="inline">-h1</code>. 
										If nothing is entered all possible nodes for the current <b>h1</b> will be considered.</td>
								</tr>
								<tr>
									<td><code class="inline">-c</code></td>
									<td>The maximum number of initial groups to consider for any gene tree. Default: 8, Max value: 18</td>
								</tr>
								<tr>
									<td><code class="inline">-o</code></td>
									<td>Output directory name. If the directory is not present, GRAMPA will created it for you.</td>
								</tr>
								<tr>
									<td><code class="inline">-f</code></td>
									<td>By default, all output files created by GRAMPA will have the prefix '<b>grampa-</b>'. You can specify a different prefix 
									with this option.</td>
								</tr>
								<tr>
									<td><code class="inline">-v</code></td>
									<td>Control the amount of output printed to the screen. 0: print nothing. 1: print only some info at the start. 2: print all log info to screen. 
									3 (default): print all log info to the screen as well as progress updates for certain steps.</td>
								</tr>
								<tr>
									<td><code class="inline">-p</code></td>
									<td>The number of processes GRAMPA should use for reconciliations.</td>
								</tr>
								<tr>
									<td><code class="inline">--multree</code></td>
									<td>Set this flag if your input species tree is a MUL-tree.</td>
								</tr>
								<tr>
									<td><code class="inline">--labeltree</code></td>
									<td>The program will simply label your input species tree.</td>
								</tr>
								<tr>
									<td><code class="inline">--numtrees</code></td>
									<td>The program will simply count the number of possible MUL-trees given <code class="inline">-s</code>. <code class="inline">-h1</code> and 
									<code class="inline">-h2</code> may also be supplied.</td>
								</tr>
								<tr>
									<td><code class="inline">--buildmultrees</code></td>
									<td>Build MUL-trees given <code class="inline">-s</code>, <code class="inline">-h1</code>, and <code class="inline">-h2</code> and write them 
									to the log file.</td>
								</tr>
								<tr>
									<td><code class="inline">--checknums</code></td>
									<td>If this flag is entered, the program will just calculate the number of groups per gene tree and exit. No reconciliations will be done.</td>
								</tr>
								<tr>
									<td><code class="inline">--st-only</code></td>
									<td>Only do reconciliations to the input singly-labeled species tree.</td>
								</tr>
								<tr>
									<td><code class="inline">--no-st</code></td>
									<td>Skip doing reconciliations to the input singly-labled species tree.</td>
								</tr>
								<tr>
									<td><code class="inline">--maps</code></td>
									<td>Output the node maps for each reconciliation in addition to the scores. The maps will be placed in the detailed output file.</td>
								</tr>
								<tr>
									<td><code class="inline">--version</code></td>
									<td>Print out version info and exit.</td>
								</tr>
								<tr>
									<td><code class="inline">--tests</code></td>
									<td>Run the tests script.</td>
								</tr>
							</tbody>	
						</table>
					</div>   

				</br>
				<div class="row">
					<a name="detail_options"></a>
					<div class="col-1" id="jump_row">
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

					<a name="-s"></a><h3><code class="inline">-s</code> : A rooted, Newick formatted species tree. This tree can be singly-labeled or MUL.</h3>
						<ul>
							<p>The tree can be in a file, in which case you enter the file name here, or you can simply paste the tree string into the command line.</p>

							<p>Entering a singly-labeled tree means you wish to search for the most parsimonious polyploidy scenario. GRAMPA will build MUL-trees based
								on this singly-labeled tree and calculate reconciliation scores. You can specify the range of MUL-trees to build with the <code class="inline">
								-h1</code> and <code class="inline">-h2</code> options.</p>

							<p>Example singly-labeled species tree:</p>

							<center><code class="cb">(((a,(x,(y,z))),b),(c,d))</code></center>

							<p>Entering a MUL-tree is the equivalent of entering a singly-labeled tree and specifying a single H1 and single H2 node. It represents a 
								single scenario of polyploidy and can be used if you wish to count the number of duplications and losses on gene trees given that
								scenario.</p>

							<p><b>NOTE: If a MUL-tree is entered, the <code class="cb">--multree</code> flag must be set.</b></p>

							<p>Example MUL-tree:</p>

							<center><code class="cb">((((a,(x,(y,z))),b),(x,(y,z))),(c,d))</code></center>
						</ul>

					<a name="-g"></a><h3><code class="inline">-g</code> : A file containing newick formatted gene trees.</h3>

						<ul>
							<p>This file should contain one or more bifurcating, Newick formatted gene trees, with one tree per line in the file. Currentky, gene trees with
								unresolved nodes (polytomies) are not supported as they falsely increase the number of losses counted in that tree.</p>
								
							<p><b>The tip labels in the gene trees must end with _[species label]</b> where [species label] matches a tip label in the species tree
								This is necessary so GRAMPA can initialize the mappings correctly.</p>

							<p>Alternatively, if you wish to reconcile to only a single gene tree, you can simply paste the tree string into the command line.</p>
						</ul>

					<a name="h1h2"></a><h3><code class="inline">-h1</code> and <code class="inline">-h2</code> : GRAMPA's search parameters.</h3>
						<ul>
							<p>H1 and H2 are nodes in the singly-labeled species tree that define how to build a MUL-tree. H1 is the node that represents the
								polyploid clade. The subtree rooted at H1 and the branch that H1 subtends will be copied onto the branch that H2 subtends:</p>

							<div class="row img-row">
								<div class="col-7-24 img-margin-left"></div>
								<div class="col-10-24 img-col">
									<img class="grid-img" src="img/mul_tree.png">
								</div>
								<div class="col-7-24 img-margin-right"></div>
							</div>

							<p>In the above example, H1 is node 2 and H2 is node 5 in the singly-labeled tree. This leads to the MUL-tree on the right.</p>
							<p>H1 and H2 can be input in 2 different, equivalent ways:</p>

							<center><code class="cb">-h1 "2" -h2 "5"</code> and <code class="cb">-h1 "x,y,z" -h2 "c,d"</code></center>

							<p>The first way relies on internal node labels. To label your species tree, use the <code class="inline">--labeltree</code> option.
								</br><b>IMPORTANT: For now, only use node labels as specified by <code class="inline">--labeltree</code>. Custom labels will not work.</b></p>

							<p>The second way uses a list of the species that define that node. <b>Species must be comma delimited.</b></p>

							<p>H2 cannot be located below H1 in the species tree! If this occurs, GRAMPA will just tell you that it's not possible and move on.</p>

							<p>Multiple H1 and H2 nodes can be entered as a space delimited list:</p>

							<center><code class="cb">-h1 "2 3" -h2 "5 6"</code> and <code class="cb">-h1 "x,y,z a,x,y,z" -h2 "c,d a,b,c,d,x,y,z"</code> are equivalent.</center>

							<p>Entering this means that GRAMPA will first set H1 as node 2 and try both nodes 5 and 6 as H2. Then H1 will be set to node 3 and will try nodes 5 
								and 6 as H2.</p>

							<h4>If <code class="inline">-h1</code> and <code class="inline">-h2</code> are not specified, GRAMPA will try all possible node combinations
								of H1 and H2!</h4>
						</ul>

					<a name="-c"></a><h3><code class="inline">-c</code> : The group cap</h3>
						<ul>
							<p>GRAMPA uses the standard LCA reconciliation algorithm on MUL-trees, meaning that some genes have more than one possible mapping.
								We get around this by trying ALL possible initial mappings and picking the one with the lowest score. This works, but also means
								our program has an exponential runtime based on the number of genes from polyploid species in any given gene tree. We get around
								this in several ways by collapsing and fixing groups (see paper), but there can still be lots of groups. This parameter sets the
								maximum number of groups to consider for any gene tree. If a gene tree has more than this number of groups, it will be skipped.</p>

							<p>Default is 8 groups, with a max setting of 18.</p>
						</ul>

					<a name="-o"></a><h3><code class="inline">-o</code> : Output directory</h3>
						<ul>
							<p>Grampa creates several output files, so it is easiest just to place them all in a single directory. That directory can be specified
								with this option, and will be created for you if it doesn't exist. If this option is not specified, the default output directory
								is "grampa_[date]-[time]".</p>
						</ul>

					<a name="-f"></a><h3><code class="inline">-f</code> : Output file prefix</h3>
						<ul>
							<p>By default, all output files created by GRAMPA will have the prefix '<b>grampa-</b>'. You can specify a different prefix with this option.
								For example, a run with <code class="inline">-f test</code> will generate the following output files, all within the output directory:</p>

							<center><code class="cb">test_out.txt, test_det.txt, test_checknums.txt</code></center>
						</ul>

					<a name="--multree"></a><h3><code class="inline">--multree</code> : Input MUL-tree flag</h3>
						<ul>
							<p>GRAMPA can accept both singly-labeled and MUL-trees as input. If your input species tree (<code class="cb">-s</code>) is a MUL-tree,
								you must set this flag so GRAMPA knows to read it as a MUL-tree. A MUL-tree represents a single possible polyploid scenario and it is
								equivalent to entering a singly-labeled tree with a single H1 and H2 node specified.</p>
						</ul>

					<a name="--labeltree"></a><h3><code class="inline">--labeltree</code> : Species tree labeling</h3>
						<ul>
							<p>This option can be used in conjunction with <code class="inline">-s</code> to simply add internal node labels to a species tree and print
								it to the screen. For example, if the file <b>species.tree</b> contains the following tree:</p>

							<center><code class="cb">(((a,(x,(y,z))),b),(c,d))</code></center>

							<p>Then the command:</p>

							<center><code class="cb">python grampa.py -s species.tree --labeltree</code></center>

							<p>Will simply print this to the screen as output:</p>

							<center><code class="cb">(((a,(x,(y,z)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(c,d)&lt;5&gt;)&lt;6&gt;</code></center>
						</ul>

					<a name="--numtrees"></a><h3><code class="inline">--numtrees</code> : Building MUL-trees</h3>
						<ul>
							<p>
								This option quickly calculates how many MUL-trees are to be built with a given H1 and H2 set. If neither H1 or H2 are set
                                it will display the total number of MUL-trees possible for the input species tree. This information is printed to the screen.
                            </p>
						</ul>                        
                        
					<a name="--buildmultrees"></a><h3><code class="inline">--buildmultrees</code> : Building MUL-trees</h3>
						<ul>
							<p>This option can be used with <code class="inline">-s</code>, <code class="inline">-h1</code>, and <code class="inline">-h2</code> to build
								MUL-trees from a standard species tree. For example, if the file <b>species.tree</b> contains the following tree:</p>

							<center><code class="cb">(((a,(x,(y,z))),b),(c,d))</code></center>

							<p>Then the command:</p>

							<center><code class="cb">python grampa.py -s species.tree -h1 "2" -h2 "4" -o multree_ex --buildmultrees</code></center>

							<p>Will yield the following output in the main output file (<b>multree_ex/grampa-out.txt</b>):</p>

							<center><code class="cb">((((a,(x+,(y+,z+)&lt;1&gt;)&lt;2&gt;)&lt;3&gt;,b)&lt;4&gt;,(x*,(y*,z*)&lt;5&gt;)&lt;6&gt;)&lt;7&gt;,(c,d)&lt;8&gt;)&lt;9&gt;</code></center>
                            
                            <p>
								The MUL-trees are written to the log file.
                            </p>
						</ul>

					<a name="--checknums"></a><h3><code class="inline">--checknums</code> : Group counting</h3>
						<ul>
							<p>With this set, the program will run normally with the specified options, except no reconciliations will be done. Instead, only the
								checknums output file will be created and will contain the number of polyploid groups for each gene tree. Use this to decide the
								best setting for <code class="inline">-c</code>.</p>
						</ul>

					<a name="--st-only"></a><h3><code class="inline">--st-only</code> : Group counting</h3>
						<ul>
							<p>
                            	By default, GRAMPA reconciles the gene trees to all specified MUL-trees as well as the singly-labeled input species tree. Set this option to ONLY
                                do reconciliations to the singly-labeled input species tree.
                            </p>
						</ul>                        

					<a name="--no-st"></a><h3><code class="inline">--no-st</code> : Group counting</h3>
						<ul>
							<p>
                            	By default, GRAMPA reconciles the gene trees to all specified MUL-trees as well as the singly-labeled input species tree. Set this option to SKIP
                                reconciliations to the singly-labeled input species tree.
                            </p>
						</ul>  
                                                
					<a name="--maps"></a><h3><code class="inline">--maps</code> : Output node mappings</h3>
						<ul>
							<p>
                            	This option adds a column to the <b>grampa-detailed.txt</b> with the actual LCA node mappings for each gene tree and MUL-tree combo. The column
                                contains a Newick formatted version of the gene tree with nodes labeled as follows:
                            </p>
                            
                            <center><code class="cb">Node[Map-Dups]</code></center>
                            
                            <p>
								Where Map indicates the node in the MUL-tree that this gene tree node maps to and Dups the number of duplications this mapping incurs. These trees 
                                can be rendered with a tree viewer such as <a href="http://doua.prabi.fr/software/seaview" target="_blank">SeaView</a>
								or <a href="http://tree.bio.ed.ac.uk/software/figtree/" target="_blank">FigTree</a>.
                            </p>
<!--
<pre><code>mul.tree        gene.tree       dups    losses  total.score     maps
106     1       1       1       2       (((1_a[a+-0],((1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;3&gt;[&lt;1&gt;+-0])&lt;4&gt;[&lt;2&gt;+-0])&lt;5&gt;[&lt;2&gt;+-1])&lt;6&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;7&gt;[&lt;4&gt;+-0],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
106     2       1       2       3       ((((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;4&gt;[&lt;1&gt;+-0])&lt;5&gt;[&lt;2&gt;+-0])&lt;6&gt;[&lt;3&gt;+-1],1_b[b+-0])&lt;7&gt;[&lt;4&gt;+-0],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
106     3       1       3       4       ((((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;4&gt;[&lt;4&gt;+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;5&gt;[&lt;1&gt;+-0])&lt;6&gt;[&lt;2&gt;+-0])&lt;7&gt;[&lt;4&gt;+-1],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
106     4       1       2       3       ((((1_a[a+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;4&gt;[&lt;1&gt;+-0])&lt;5&gt;[&lt;2&gt;+-0])&lt;6&gt;[&lt;3&gt;+-1],1_b[b+-0])&lt;7&gt;[&lt;4&gt;+-0],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
106     5       1       3       4       (((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],(1_b[b+-0],(2_x[x+-0],(2_y[y+-0],2_z[z+-0])&lt;4&gt;[&lt;1&gt;+-0])&lt;5&gt;[&lt;2&gt;+-0])&lt;6&gt;[&lt;4&gt;+-0])&lt;7&gt;[&lt;4&gt;+-1],(1_c[c+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
106     6       1       3       4       (((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;4&gt;[&lt;4&gt;+-0],((1_c[c+-0],(2_x[x*-0],(2_y[y*-0],2_z[z*-0])&lt;5&gt;[&lt;5&gt;+-0])&lt;6&gt;[&lt;6&gt;+-0])&lt;7&gt;[&lt;8&gt;+-0],1_d[d+-0])&lt;8&gt;[&lt;8&gt;+-1])&lt;9&gt;[&lt;9&gt;+-0]
106     7       0       0       0       (((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;4&gt;[&lt;4&gt;+-0],(1_c[c+-0],(1_d[d+-0],(2_x[x*-0],(2_y[y*-0],2_z[z*-0])&lt;5&gt;[&lt;5&gt;+-0])&lt;6&gt;[&lt;6&gt;+-0])&lt;7&gt;[&lt;7&gt;+-0])&lt;8&gt;[&lt;8&gt;+-0])&lt;9&gt;[&lt;9&gt;+-0]
106     8       1       3       4       (((1_a[a+-0],(1_x[x+-0],(1_y[y+-0],1_z[z+-0])&lt;1&gt;[&lt;1&gt;+-0])&lt;2&gt;[&lt;2&gt;+-0])&lt;3&gt;[&lt;3&gt;+-0],1_b[b+-0])&lt;4&gt;[&lt;4&gt;+-0],((1_c[c+-0],1_d[d+-0])&lt;5&gt;[&lt;8&gt;+-0],(2_x[x*-0],(2_y[y*-0],2_z[z*-0])&lt;6&gt;[&lt;5&gt;+-0])&lt;7&gt;[&lt;6&gt;+-0])&lt;8&gt;[&lt;8&gt;+-1])&lt;9&gt;[&lt;9&gt;+-0]
</code></pre>
-->							
						</ul>
			</div>
		</div>
		<div class="col-3-24" id="margin"></div>
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
nav = RC.readNav(pagefile, "", "results/yeast/", "results/wheat/");
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));