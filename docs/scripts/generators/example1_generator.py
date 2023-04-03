############################################################
# For GRAMPA site, 12.19
# This generates the file "example1.html"
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

				<h1>Infer presence and mode of polyploidy</h1>

				<h4>Below are the inputs, commands, and outputs to do several analyses with GRAMPA. The inputs are based on
					simulated data. For more detailed info on the simulations check
					<a href="http://biorxiv.org/content/early/2017/03/21/058149" target="_blank">our paper</a>.</h4>

				<div id="msg-cont">
					<div id="msg">
						<div id="msg-banner">Important!</div>
						<div id="msg-text">
							<p>
								These examples were done with earlier versions of GRAMPA (<1.4.0) so some of the command line options and output formats
								may have changed, but the general idea and results remain the same. See the <a href="readme.html">README</a> for up-to-date info
								on options and formats.
							</p>
							<p></p>
						</div>
					</div>
				</div>
				<div class="sep-div-2"></div>  

				<div class="row">
					<a name="allo"></a>
					<div class="col-1" id="jump_row">
						<center>
						<div id="jump_container">
							Jump to section:
							<a class="jump_link" href="#allo">Allopolyploidy</a>
							<a class="jump_link" href="#auto">Autopolyploidy</a>
							<a class="jump_link" href="#nop">No polyploidy</a>
						</div>
						</center>
					</div>
				</div>	

				<h1>Allopolyploidy</h1>

				<h3>Inputs</h3>

				<p>1000 gene trees were simulated with gain and loss using <a href="https://github.com/arvestad/jprime" target="_blank">JPrIME</a> based on the
					following allopolyploid-like MUL-tree:</p>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/allo_3a/allo_m.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>				

				<p>In this scenario, lineages B and C hybridized to form an allopolyploid lineage that diversified into the x,y,z clade.</p>

				<p>We then remove one polyploid clade from the MUL-tree to get a singly-labeled tree as input for GRAMPA. This is the same topology as above, except
					that the x,y,z clade sister to C is removed. This is the type of tree that typical phylogenetic reconstruction programs would produce even in
					the presence of allopolyploidy:</p>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/allo_3a/allo_s.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>				

				<p>So, the input files for this search are:</p>

				<ol>
					<li>Singly-labeled species tree: <a href="example_data/ex1/allo_3a/spec_tree_3a.tre" download>spec_tree_3a.tre</a></li>
					<li>1000 gene trees simulated from an allopolyploid MUL-tree: <a href="example_data/ex1/allo_3a/gene_trees_3a.txt" download>gene_trees_3a.txt</a></li>
				</ol>

				<h3>GRAMPA command</h3>

				<p>Since in reality we wouldn't know whether there is an allo-, auto-, or no polyploidy in this tree, we want GRAMPA to search all nodes as possible
					polyploid lineages. That means we don't specify <code class="inline">-h1</code> or <code class="inline">-h2</code>.

				<center><pre class="cmd"><code>python grampa.py -s spec_tree_3a.tre -g gene_trees_3a.txt -o allo_example_output -f allo_test -v 0</code></pre></center>

				<h3>Outputs</h3>

				<p>The above command would create the directory <code class="cb">allo_example_output</code> with four output files</p>

				<ul>
					<li><a href="example_data/ex1/allo_3a/allo_example_output/allo_test_checknums.txt" download>allo_test_checknums.txt</a></li>
					<li><a href="example_data/ex1/allo_3a/allo_example_output/allo_test_det.txt" download>allo_test_det.txt</a></li>
					<li><a href="example_data/ex1/allo_3a/allo_example_output/allo_test_out.txt" download>allo_test_out.txt</a></li>
					<li><a href="example_data/ex1/allo_3a/allo_example_output/allo_test_trees_filtered.txt" download>allo_test_trees_filtered.txt</a></li>
				</ul>

				<p>Since we are trying to determine the mode of polyploidy, we are interested in the <code class="cb">allo_test_out.txt</code> file. This file contains
					log info and the total reconciliation scores for each MUL-tree considered and looks something like this:</p>

<pre><code># Tree #    H1 node H2 node Tree string Total score
ST          ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt; 7980
MT-1    A   A   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,(A+,A*)&lt;4&gt;)&lt;5&gt;,C)&lt;6&gt;,D)&lt;7&gt;    8272
MT-2    A   C   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A+)&lt;4&gt;,(C,A*)&lt;5&gt;)&lt;6&gt;,D)&lt;7&gt;    8767
MT-3    A   B   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,(B,A*)&lt;3&gt;)&lt;4&gt;,A+)&lt;5&gt;,C)&lt;6&gt;,D)&lt;7&gt;    8777
MT-4    A   D   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A+)&lt;4&gt;,C)&lt;5&gt;,(D,A*)&lt;6&gt;)&lt;7&gt;    8553
.
.
.
MT-126  &lt;5&gt; &lt;6&gt; (((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B+)&lt;3&gt;,A+)&lt;4&gt;,C+)&lt;5&gt;,D)&lt;6&gt;,(((((x*,y*)&lt;7&gt;,z*)&lt;8&gt;,B*)&lt;9&gt;,A*)&lt;10&gt;,C*)&lt;11&gt;)&lt;12&gt;    8420
MT-127  &lt;5&gt; &lt;5&gt; (((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B+)&lt;3&gt;,A+)&lt;4&gt;,C+)&lt;5&gt;,(((((x*,y*)&lt;6&gt;,z*)&lt;7&gt;,B*)&lt;8&gt;,A*)&lt;9&gt;,C*)&lt;10&gt;)&lt;11&gt;,D)&lt;12&gt;    7824
# ---------
The MUL-tree with the minimum parsimony score is MT-74: ((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,(C,((x*,y*)&lt;5&gt;,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;
Score = 5018
</code></pre>

				<p>GRAMPA tells us MUL-tree 74 is the lowest scoring tree:</p>

				<center><code class="cb">((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,(C,((x*,y*)&lt;5&gt;,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;</code></center>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/allo_3a/mt_74.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<p>Notice that this is the same topology that was used to simulate the gene-trees. GRAMPA has successfully identified an allopolyploid MUL-tree and
					placed the second polyploid lineage on the correct branch!</p>

				<div class="row">
					<a name="auto"></a>
					<div class="col-1" id="jump_row">
						<center>
						<div id="jump_container">
							Jump to section:
							<a class="jump_link" href="#allo">Allopolyploidy</a>
							<a class="jump_link" href="#auto">Autopolyploidy</a>
							<a class="jump_link" href="#nop">No polyploidy</a>
						</div>
						</center>
					</div>
				</div>	

				<h1>Autopolyploidy</h1>

				<p>1000 gene trees were simulated with gain and loss using <a href="https://github.com/arvestad/jprime" target="_blank">JPrIME</a> based on the
					following autopolyploid-like MUL-tree:</p>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/auto_18/auto_m.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<p>In this scenario, a lineage sister to species C underwent autopolyploidization and subsequently diversified into the x,y,z clade.</p>

				<p>We then remove one polyploid clade from the MUL-tree to get a singly-labeled tree as input for GRAMPA. This is the same topology as above, except
					that one x,y,z clade is removed. This is the type of tree that typical phylogenetic reconstruction programs would produce even in
					the presence of autopolyploidy:</p>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/auto_18/auto_s.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<p>So, the input files for this search are:</p>

				<ol>
					<li>Singly-labeled species tree: <a href="example_data/ex1/auto_18/spec_tree_18.tre" download>spec_tree_18.tre</a></li>
					<li>1000 gene trees simulated from an allopolyploid MUL-tree: <a href="example_data/ex1/auto_18/gene_trees_18.txt" download>gene_trees_18.txt</a></li>
				</ol>

				<h3>GRAMPA command</h3>

				<p>Since in reality we wouldn't know whether there is an allo-, auto-, or no polyploidy in this tree, we want GRAMPA to search all nodes as possible
					polyploid lineages. That means we don't specify <code class="cb">-h1</code> or <code class="cb">-h2</code>.

				<center><pre class="cmd"><code>python grampa.py -s spec_tree_18.tre -g gene_trees_18.txt -o auto_example_output -f auto_test -v 0</code></pre></center>

				<h3>Outputs</h3>

				<p>The above command would create the directory <code class="cb">auto_example_output</code> with four output files</p>

				<ul>
					<li><a href="example_data/ex1/auto_18/auto_example_output/auto_test_checknums.txt" download>auto_test_checknums.txt</a></li>
					<li><a href="example_data/ex1/auto_18/auto_example_output/auto_test_det.txt" download>auto_test_det.txt</a></li>
					<li><a href="example_data/ex1/auto_18/auto_example_output/auto_test_out.txt" download>auto_test_out.txt</a></li>
					<li><a href="example_data/ex1/auto_18/auto_example_output/auto_test_trees_filtered.txt" download>auto_test_trees_filtered.txt</a></li>
				</ul>

				<p>Since we are trying to determine the mode of polyploidy, we are interested in the <code class="cb">auto_test_out.txt</code> file. This file contains
					log info and the total reconciliation scores for each MUL-tree considered and looks something like this:</p>

<pre><code># Tree #    H1 node H2 node Tree string Total score
ST          (((B,A)&lt;1&gt;,(((x,y)&lt;2&gt;,z)&lt;3&gt;,C)&lt;4&gt;)&lt;5&gt;,D)&lt;6&gt; 5476
MT-1    A   A   (((B,(A+,A*)&lt;1&gt;)&lt;2&gt;,(((x,y)&lt;3&gt;,z)&lt;4&gt;,C)&lt;5&gt;)&lt;6&gt;,D)&lt;7&gt;    6280
MT-2    A   C   (((B,A+)&lt;1&gt;,(((x,y)&lt;2&gt;,z)&lt;3&gt;,(C,A*)&lt;4&gt;)&lt;5&gt;)&lt;6&gt;,D)&lt;7&gt;    6244
MT-3    A   B   ((((B,A*)&lt;1&gt;,A+)&lt;2&gt;,(((x,y)&lt;3&gt;,z)&lt;4&gt;,C)&lt;5&gt;)&lt;6&gt;,D)&lt;7&gt;    6115
MT-4    A   D   (((B,A+)&lt;1&gt;,(((x,y)&lt;2&gt;,z)&lt;3&gt;,C)&lt;4&gt;)&lt;5&gt;,(D,A*)&lt;6&gt;)&lt;7&gt;    6088
.
.
.
MT-132  &lt;5&gt; &lt;6&gt; ((((B+,A+)&lt;1&gt;,(((x+,y+)&lt;2&gt;,z+)&lt;3&gt;,C+)&lt;4&gt;)&lt;5&gt;,D)&lt;6&gt;,((B*,A*)&lt;7&gt;,(((x*,y*)&lt;8&gt;,z*)&lt;9&gt;,C*)&lt;10&gt;)&lt;11&gt;)&lt;12&gt;    6040
MT-133  &lt;5&gt; &lt;5&gt; ((((B+,A+)&lt;1&gt;,(((x+,y+)&lt;2&gt;,z+)&lt;3&gt;,C+)&lt;4&gt;)&lt;5&gt;,((B*,A*)&lt;6&gt;,(((x*,y*)&lt;7&gt;,z*)&lt;8&gt;,C*)&lt;9&gt;)&lt;10&gt;)&lt;11&gt;,D)&lt;12&gt;    6103
# ---------
The MUL-tree with the minimum parsimony score is MT-57: (((B,A)&lt;1&gt;,((((x+,y+)&lt;2&gt;,z+)&lt;3&gt;,((x*,y*)&lt;4&gt;,z*)&lt;5&gt;)&lt;6&gt;,C)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;
Score = 4807
</code></pre>

				<p>GRAMPA tells us MUL-tree 57 is the lowest scoring tree:</p>

				<center><code class="cb">(((B,A)&lt;1&gt;,((((x+,y+)&lt;2&gt;,z+)&lt;3&gt;,((x*,y*)&lt;4&gt;,z*)&lt;5&gt;)&lt;6&gt;,C)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;</code></center>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/auto_18/mt_57.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<p>Notice that this is the same topology that was used to simulate the gene-trees. GRAMPA has successfully identified an autopolyploid MUL-tree on the
					correct branch!</p>

				<div class="row">
					<a name="nop"></a>
					<div class="col-1" id="jump_row">
						<center>
						<div id="jump_container">
							Jump to section:
							<a class="jump_link" href="#allo">Allopolyploidy</a>
							<a class="jump_link" href="#auto">Autopolyploidy</a>
							<a class="jump_link" href="#nop">No polyploidy</a>
						</div>
						</center>
					</div>
				</div>	

				<h1>No polyploidy</h1>

				<p>1000 gene trees were simulated with gain and loss <a href="https://github.com/arvestad/jprime" target="_blank">JPrIME</a> based on the
					following singly-labeled tree:</p>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/nop_33/nop_s.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<p>In this scenario, no polyploidy has occurred and this is the same tree we give to GRAMPA.</p>

				<p>So, the input files for this search are:</p>

				<ol>
					<li>Singly-labeled species tree: <a href="example_data/ex1/nop_33/spec_tree_33.tre" download>spec_tree_33.tre</a></li>
					<li>1000 gene trees simulated from an allopolyploid MUL-tree: <a href="example_data/ex1/nop_33/gene_trees_33.txt" download>gene_trees_33.txt</a></li>
				</ol>

				<h3>GRAMPA command</h3>

				<p>Since in reality we wouldn't know whether there is an allo-, auto-, or no polyploidy in this tree, we want GRAMPA to search all nodes as possible
					polyploid lineages. That means we don't specify <code class="inline">-h1</code> or <code class="inline">-h2</code>.

				<center><pre class="cmd"><code>python grampa.py -s spec_tree_33.tre -g gene_trees_33.txt -o nop_example_output -f nop_test -v 0</code></pre></center>

				<h3>Outputs</h3>

				<p>The above command would create the directory <code class="cb">nop_example_output</code> with four output files</p>

				<ul>
					<li><a href="example_data/ex1/nop_33/nop_example_output/nop_test_checknums.txt" download>nop_test_checknums.txt</a></li>
					<li><a href="example_data/ex1/nop_33/nop_example_output/nop_test_det.txt" download>nop_test_det.txt</a></li>
					<li><a href="example_data/ex1/nop_33/nop_example_output/nop_test_out.txt" download>nop_test_out.txt</a></li>
					<li><a href="example_data/ex1/nop_33/nop_example_output/nop_test_trees_filtered.txt" download>nop_test_trees_filtered.txt</a></li>
				</ul>

				<p>Since we are trying to determine the mode of polyploidy, we are interested in the <code class="cb">nop_test_out.txt</code> file. This file contains
					log info and the total reconciliation scores for each MUL-tree considered and looks something like this:</p>

<pre><code># Tree #    H1 node H2 node Tree string Total score
ST          ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt; 4115
MT-1    A   A   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,(A+,A*)&lt;4&gt;)&lt;5&gt;,C)&lt;6&gt;,D)&lt;7&gt;    4423
MT-2    A   C   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A+)&lt;4&gt;,(C,A*)&lt;5&gt;)&lt;6&gt;,D)&lt;7&gt;    4753
MT-3    A   B   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,(B,A*)&lt;3&gt;)&lt;4&gt;,A+)&lt;5&gt;,C)&lt;6&gt;,D)&lt;7&gt;    4929
MT-4    A   D   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A+)&lt;4&gt;,C)&lt;5&gt;,(D,A*)&lt;6&gt;)&lt;7&gt;    4754
.
.
.
MT-126  &lt;5&gt; &lt;6&gt; (((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B+)&lt;3&gt;,A+)&lt;4&gt;,C+)&lt;5&gt;,D)&lt;6&gt;,(((((x*,y*)&lt;7&gt;,z*)&lt;8&gt;,B*)&lt;9&gt;,A*)&lt;10&gt;,C*)&lt;11&gt;)&lt;12&gt;    4760
MT-127  &lt;5&gt; &lt;5&gt; (((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B+)&lt;3&gt;,A+)&lt;4&gt;,C+)&lt;5&gt;,(((((x*,y*)&lt;6&gt;,z*)&lt;7&gt;,B*)&lt;8&gt;,A*)&lt;9&gt;,C*)&lt;10&gt;)&lt;11&gt;,D)&lt;12&gt;    4877
# ---------
The tree with the minimum parsimony score is the singly-labled tree (ST):   ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt;
Score = 4115
</code></pre>

				<p>GRAMPA tells us that the singly-labeled tree is the lowest scoring tree:</p>

				<center><code class="cb">((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt;</code></center>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex1/nop_33/nop_s.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<p>Notice that this is the same topology that was used to simulate the gene-trees. GRAMPA has successfully determined that no polyploidy has
					occurred among these lineages!</p>
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
pagefile = "example1.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Identify"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile, "", "results/yeast/", "results/wheat/");
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));