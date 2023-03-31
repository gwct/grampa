############################################################
# For GRAMPA site, 12.19
# This generates the file "example2.html"
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

				<h1>Placing a known WGD on a phylogeny</h1>

				<h4>Below are the inputs, commands, and outputs to do an analysis with GRAMPA to place a known WGD on a phylogeny. The inputs are based on
					simulated data. For more detailed info on the simulations check
					<a href="http://biorxiv.org/content/early/2017/03/21/058149" target="_blank">our paper</a>.</h4>

				<h3>Inputs</h3>

				<p>Suppose you have a set of species, of which you have evidence that a some may be the result of a polyploidization event. You also may have
					an idea about the parental lineages of the polyploid species. So you build a species tree of your taxa and, since species tree reconstruction
					programs output singly-labeled trees, you get this as the result:</p>

				<div class="row img-row">
					<div class="col-8-24 img-margin-left"></div>
					<div class="col-8-24 img-col">
						<img class="grid-img" src="example_data/ex2/spec_tree.png">
					</div>
					<div class="col-8-24 img-margin-right"></div>
				</div>			

				<p>Your hypothesis is that species x, y, and z may be the result of polyploidization (but you're not sure if all of them are). The singly-labeled tree
					implicitly identifies one of the parental lineages of the suspected polyploids by placing the polyploid species sister to it (in this case lineage
					B seems to be one of the parents). You also think this tree may be the result of an allopolyploidy, and you think some lineage sister to C, D, the
					C,D clade, or at the root of the tree may have hybridized with some lineage related to B to form the x,y,z clade. Working with this prior knowledge
					GRAMPA can test your hypotheses of polyploidization.</p>

				<p>The input files you would need are:</p>

				<ol>
					<li>Singly-labeled species tree: <a href="example_data/ex2/spec_tree_3a.tre" download>spec_tree_3a.tre</a></li>
					<li>Gene trees from your set of species (in this case 1000 gene trees simulated with gain and loss) : 
						<a href="example_data/ex2/gene_trees_3a.txt" download>gene_trees_3a.txt</a></li>
				</ol>

				<h3>GRAMPA command</h3>

				<p>Since we have some idea of the lineages involved in the polyploidization event, we would want to limit GRAMPA's search to those lineages
					with the <code class="cb">-h1</code> and <code class="cb">-h2</code> search parameters.</p>

				<center><pre class="cmd"><code>python grampa.py -s spec_tree_3a.tre -g gene_trees_3a.txt -h1 "x 1 2" -h2 "C D 5 6" -o ex2_output -f ex2_test -v 0</code></pre></center>
				
				<p>Above we have specified <code class="inline">-h1</code> and <code class="inline">-h2</code> by using the node labels in the tree. Alternatively, we could specify
					an equivalent <code class="inline">-h1</code> and <code class="inline">-h2</code> search by defining the labels based on the sets of tips that define them:</p>

				<center><pre class="cmd"><code >python grampa.py -s spec_tree_3a.tre -g gene_trees_3a.txt -h1 "x x,y x,y,z" -h2 "C D A,x,y,z,B,C A,x,y,z,B,C,D" -o ex2_output -f ex2_test -v 0</code></pre></center>
				
				<p>These two commands are equivalent. The second method is slightly more cumbersome, but does not require you to have internal labels on your tree.
					Although, GRAMPA can easily add internal labels to your input tree with the <code class="inline">--labeltree</code> command.</p>

				<h3>Outputs</h3>

				<p>The above command would create the directory <b>ex2_output</b> with three output files</p>
				<ul>
					<li><a href="example_data/ex2/ex2_output/ex2_test_checknums.txt" download>ex2_test_checknums.txt</a></li>
					<li><a href="example_data/ex2/ex2_output/ex2_test_det.txt" download>ex2_test_det.txt</a></li>
					<li><a href="example_data/ex2/ex2_output/ex2_test_out.txt" download>ex2_test_out.txt</a></li>
				</ul>
				<p>Since we are trying to determine the mode of polyploidy, we are interested in the <b>ex2_test_out.txt</b> file. This file contains
					log info and the total reconciliation scores for each MUL-tree considered and looks something like this:</p>

<pre><code># Tree #    H1 node H2 node Tree string Total score
ST          ((((((x,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt; 8312
MT-1    x   C   ((((((x+,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,(C,x*)&lt;5&gt;)&lt;6&gt;,D)&lt;7&gt;    9109
MT-2    x   D   ((((((x+,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,(D,x*)&lt;6&gt;)&lt;7&gt;    8949
MT-3    x   &lt;5&gt; (((((((x+,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,x*)&lt;6&gt;,D)&lt;7&gt;    9259
MT-4    x   &lt;6&gt; (((((((x+,y)&lt;1&gt;,z)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt;,x*)&lt;7&gt;    9304
.
.
.
MT-11   &lt;2&gt; &lt;5&gt; (((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,((x*,y*)&lt;6&gt;,z*)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;  7598
MT-12   &lt;2&gt; &lt;6&gt; (((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,C)&lt;5&gt;,D)&lt;6&gt;,((x*,y*)&lt;7&gt;,z*)&lt;8&gt;)&lt;9&gt;  8845
# ---------
The MUL-tree with the minimum parsimony score is MT-9:  ((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,(C,((x*,y*)&lt;5&gt;,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;
Score = 5242
</code></pre>

				<p>GRAMPA tells us MUL-tree 9 is the lowest scoring tree:</p>

				<center><code class="cb">((((((x+,y+)&lt;1&gt;,z+)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,(C,((x*,y*)&lt;5&gt;,z*)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;</code></center>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex2/mt_9.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>			

				<p>In other words, GRAMPA has identified the x,y,z clade as the polyploid clade and has identified the C lineage as the second parental lineage!</p>
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
pagefile = "example2.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Place"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile, "", "results/yeast/", "results/wheat/");
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));