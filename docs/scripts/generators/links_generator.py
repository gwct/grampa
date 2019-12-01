############################################################
# For GRAMPA site, 12.19
# This generates the file "links.html"
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
				<h1>Links</h1>
				<h3>People and places related to GRAMPA:</h3>
					<ul>
						<p><a href="http://www.indiana.edu/~hahnlab/" target="_blank">Hahn lab</a> / GRAMPA was developed in Matthew Hahn's lab at Indiana University</p>
						<p><a href="https://github.com/gwct/grampa" target="_blank">GRAMPA github</a> / GRAMPA's github repository.</p>	
					</ul>
				<h3>Helpful software:</h3>
					<ul>
						<p><a href="http://www.drive5.com/muscle/" target="_blank">MUSCLE</a> / Fast alignment software for DNA or amino acid sequences.</p>
						<p><a href="http://sco.h-its.org/exelixis/web/software/raxml/index.html" target="_blank">RAxML</a> / Maximum likelihood software for making
							gene trees from DNA or amino acid alignments.</p>
						<p><a href="https://github.com/smirarab/ASTRAL" target="_blank">ASTRAL</a> / Coalescent based software for making species trees from gene trees.</p>
						<p><a href="http://www.cs.cmu.edu/~durand/Notung/" target="_blank">Notung</a> / Another reconciliation program with lots of nice functions.
							Can be used to root gene trees based on a DupLoss model and for bootstrap rearrangement.</p>
						<p><a href="http://tree.bio.ed.ac.uk/software/figtree/" target="_blank">FigTree</a> / Software for rendering phylogenetic trees. All trees
							in the Examples and Results sections were rendered in FigTree and edited with Adobe Illustrator.</p>
						<p><a href="http://doua.prabi.fr/software/seaview" target="_blank">SeaView</a> / A versatile program for quick viewing of alignments and
							phylogenetic trees</p>
						<p><a href="https://github.com/arvestad/jprime" target="_blank">JPrIME</a> / Software used to simulate gene trees with varying rates of
							gene gain and loss.</p>
					</ul>
				<h3>Previous versions of GRAMPA</h3>
				<ul>
					<p><a href="prev/v1.1.zip" download>Version 1.1</a> / March 22, 2016 / Implemented gene tree filtering and several useful options. No parallelization.</p>
					<p><a href="prev/v1.0.zip" download>Version 1.0</a> / Summer 2015 / The first release. No gene tree filtering.</p>
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
pagefile = "links.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Links"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));