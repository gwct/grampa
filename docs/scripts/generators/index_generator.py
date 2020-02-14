############################################################
# For GRAMPA site, 12.19
# This generates the file "index.html"
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
				<img class="pure-img" id="logo_main" src="img/logo2.png">
				<h1>GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis</h1>
				<p>GRAMPA is a phylogenetics program for inferring the presence and mode of polyploidy in a set of species, placing a whole genome duplication event on 
				a phylogeny, and counting gene duplications and losses in the presence of polyploidy.</p>
				<p id="paper">Thomas GWC, Ather SH, Hahn MW. 2017. Gene tree reconciliation with MUL-trees to resolve polyploid analysis. <em>Systematic Biology.</em>
					66(6):1007-1018. <a href="https://doi.org/10.1093/sysbio/syx044" target="_blank">10.1093/sysbio/syx044</a></p>
				<div id="buttons_container">
					<a class="button-secondary pure-button" id="install_btn" href="about.html">About &raquo;</a><span id="buffer"></span>
					<a class="button-secondary pure-button" id="install_btn" href="readme.html">README &raquo;</a><span id="buffer"></span>
					<a class="button-secondary pure-button" id="install_btn" href="https://github.com/gwct/grampa/releases/latest" target="_blank">Download &raquo;</a>
				</div>
			</div>
		</div>
		<div class="pure-u-3-24" id="margin"></div>
	</div>

    <!--<div class="pure-g" id="q_row">
		<div class="pure-u-1-5" id="margin"></div>
		<div class="pure-u-3-5" id="q_col"><h1>What would you like to do?</h1></div>
		<div class="pure-u-1-5" id="margin"></div>
	</div> -->

	<div class="pure-g" id="task_row">
		<div class="pure-u-3-24" id="margin"></div>
		<div class="pure-u-18-24" id="task_col">
			<div class="pure-g" id="task_row_inner">
				<div class="pure-u-1-3" id="task">
					<h2>Identify the mode of polyploidy</h2>
					<p id="ex_p">GRAMPA may be able to differentiate between the cases of allopolyploidy, autopolyploidy, and no polyploidy.</p>
				</div>
				<div class="pure-u-1-3" id="task">
					<h2>Place a WGD on a phylogeny</h2>
					<p id="ex_p">If you know a WGD has taken place on your phylogeny, GRAMPA can place it on the tree.</p>
				</div>
				<div class="pure-u-1-3" id="task">
					<h2>Count duplications and losses</h2>
					<p id="ex_p">GRAMPA can accurately count duplications and losses, regardless of the presence of polyploidy.</p>
				</div>
			</div>
			<div class="pure-g" id="btn_row_inner">
				<div class="pure-u-1-3" id="task_btn"><a class="button-secondary pure-button" href="example1.html">Example &raquo;</a></div>
				<div class="pure-u-1-3" id="task_btn"><a class="button-secondary pure-button" href="example2.html">Example &raquo;</a></div>
				<div class="pure-u-1-3" id="task_btn"><a class="button-secondary pure-button" href="example3.html">Example &raquo;</a></div>
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
pagefile = "index.html";
print("Generating " + pagefile + "...");
title = "GRAMPA"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));