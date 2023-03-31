############################################################
# For GRAMPA site, 12.19
# This generates the file "example3.html"
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

				<h1>Count duplications and losses in the presence of polyploidy</h1>

				<h4>Below are the inputs, commands, and outputs to do an analysis with GRAMPA to count the number of duplications and losses in a set of gene trees
					in the presence of polyploidy. The inputs are based on simulated data. For more detailed info on the simulations check
					<a href="http://biorxiv.org/content/early/2017/03/21/058149" target="_blank">our paper</a>.</h4>

				<h3>Inputs</h3>

				<p>Suppose a hybridization event between the B and C lineages leading to the allopolyploid x,y,z clade has been identified. The MUL-tree representing this
					scenario is:</p>

				<div class="row img-row">
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">
						<img class="grid-img" src="example_data/ex3/mul_tree.png">
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>								

				<p>GRAMPA can accurately count duplications and losses with this MUL-tree as input.</p>

				<ol>
					<li>MUL-tree: <a href="example_data/ex3/mul_tree_74_3a.tre" download>mul_tree_74_3a.tre</a></li>
					<li>Gene trees from your set of species (in this case 1000 gene trees simulated with gain and loss) : 
						<a href="example_data/ex3/gene_trees_3a.txt" download>gene_trees_3a.txt</a></li>
				</ol>

				<h3>GRAMPA command</h3>

				<p>With a known polyploidy scenario (MUL-tree), we can simply reconcile to that MUL-tree:</p>

				<center><pre class="cmd"><code>python grampa.py -s mul_tree_74_3a.tre -g gene_trees_3a.txt -o ex3_output -f count_test -v 0 --multree</code></pre></center>

				<p>The <code class="inline">--multree</code> flag is required in this case to let GRAMPA know that the input species tree is a MUL-tree.</p>

				<h3>Outputs</h3>

				<p>The above command would create the directory <b>ex3_output</b> with three output files</p>

				<ul>
					<li><a href="example_data/ex3/ex3_output/count_test_checknums.txt" download>count_test_checknums.txt</a></li>
					<li><a href="example_data/ex3/ex3_output/count_test_det.txt" download>count_test_det.txt</a></li>
					<li><a href="example_data/ex3/ex3_output/count_test_out.txt" download>count_test_out.txt</a></li>
				</ul>

				<p>Since we are trying to count duplications and losses, we are interested in the <b>count_test_det.txt</b> file. This file contains
					reconciliation scores for each gene tree to the lowest scoring MUL-tree (in this case, the only MUL-tree). The contents of the file look something
					like this:</p>

<pre><code># MT-1:((((((x*,y*)&lt;1&gt;,z*)&lt;2&gt;,B)&lt;3&gt;,A)&lt;4&gt;,(C,((x+,y+)&lt;5&gt;,z+)&lt;6&gt;)&lt;7&gt;)&lt;8&gt;,D)&lt;9&gt;
# GT/MT combo   # dups  # losses    Total score
GT-1 to MT-1    5   2   7
GT-2 to MT-1    8   4   12
GT-3 to MT-1    3   0   3
GT-4 to MT-1    1   4   5
.
.
.
GT-999 to MT-1  1   2   3
GT-1000 to MT-1 2   4   6
# Gene trees with multiple maps:    36
# Total parsimony score for MT-1: 5242
</code></pre>

				<p>The total score is reported here and in the main output file (<b>count_test_out.txt</b>), but with this output you can see the exact number of duplications and
					losses for each gene tree, while properly counting multiple copies of the polyploid species (x,y,z) as either paralogs or homoeologs. Additionally, you
					could add the <code class="inline">--maps</code> flag to the command above to add another column to the detailed output file that shows the maps and
					duplication nodes in each gene tree (see the README section <a href="readme.html#--maps">--maps</a> for more info).</p>
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
pagefile = "example3.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Count"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile, "", "results/yeast/", "results/wheat/");
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));