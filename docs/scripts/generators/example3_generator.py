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
					<a href="https://doi.org/10.1093/sysbio/syx044" target="_blank">our paper</a>.</h4>
                    
				<div id="msg-cont">
					<div id="msg">
						<div id="msg-banner">Note</div>
						<div id="msg-text">
							<p>
								The examples below call GRAMPA as <code class="inline">grampa</code> assuming it has been installed from bioconda. If you installed
								from source, see <a href="readme.html#source">usage details in the README</a>.
							</p>
							<p></p>
						</div>
					</div>
				</div>                
                
				<div class="sep-div-2"></div>  

				<h3>Inputs</h3>

				<p>Suppose a hybridization event between the B and C lineages leading to the allopolyploid x,y,z clade. The MUL-tree representing this
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

				<center><pre class="cmd"><code>grampa -s mul_tree_74_3a.tre -g gene_trees_3a.txt -o ex3-output -f count-test --multree</code></pre></center>

				<p>The <code class="inline">--multree</code> flag is required in this case to let GRAMPA know that the input species tree is a MUL-tree.</p>

				<h3>Outputs</h3>

				<p>The above command would create the directory <code class="cb">ex3-output</code> with five output files</p>

				<ul>
					<li><a href="example_data/ex3/ex-3output/count-test-checknums.txt" download>count-test-checknums.txt</a></li>
					<li><a href="example_data/ex3/ex3-output/count-test-detailed.txt" download>count-test-detailed.txt</a></li>
                    <li><a href="example_data/ex3/ex3-output/count-test-dup-counts.txt" download>count-test-dup-counts.txt</a></li>
					<li><a href="example_data/ex3/ex3-output/count-test.log" download>count-test.log</a></li>
                    <li><a href="example_data/ex3/ex3-output/count-test-scores.txt" download>count-test-scores.txt</a></li>
				</ul>

				<p>Since we are trying to count duplications and losses, we are interested in the <code class="cb">count-test-detailed.txt</code> file. This file contains
					reconciliation scores for each gene tree to the lowest scoring MUL-tree (in this case, the only MUL-tree). The contents of the file look something
					like this:</p>

<pre><code>mul.tree        gene.tree       dups    losses  total.score
1       1       5       2       7
1       2       8       4       12
1       3       3       0       3
1       4       1       4       5
</code></pre>

				<p>Here you can see the exact number of duplications and losses for each gene tree when mapped to the lowest scoring MUL-tree, while properly 
                	counting multiple copies of the polyploid species (x,y,z) as either paralogs or homoeologs. Additionally, you could add the 
                    <code class="inline">--maps</code> flag to the command above to add another column to the detailed output file that shows the maps and
					duplication nodes in each gene tree (see the README section <a href="readme.html#--maps">--maps</a> for more info).
                </p>
                
                <p>Additionally, GRAMPA counts the total number of duplications across all gene trees for each node in the six lowest scoring trees. In this
                 	example, we only searched one tree so it only provides the duplication counts for that tree. These counts can
					be found in the <code class="cb">count-test-dup-counts.txt</code> file:
                </p>
                
<pre><code>mul.tree        node    dups
1       x*      102
1       y*      88
1       z*      193
1       B       234
1       A       344
1       C       232
1       x+      104
1       y+      87
1       z+      212
1       D       386
1       <1>     94
1       <2>     82
1       <3>     102
1       <4>     73
1       <5>     110
1       <6>     101
1       <7>     143
1       <8>     108
1       <9>     0
</code></pre>                
                
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