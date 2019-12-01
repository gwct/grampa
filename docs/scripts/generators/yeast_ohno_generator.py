############################################################
# For GRAMPA site, 12.19
# This generates the file "yeast_ohno.html"
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

					<li class="pure-menu-item"><a href="trees/yeast_spec.tre" class="nav_link" download>Species tree</a></li>
					<li class="pure-menu-item"><a href="trees/yeast_ohno_trees.txt" class="nav_link" download>Gene trees</a></li>
					<li class="pure-menu-item"><a href="output/yeast_ohno_out.txt" class="nav_link" download>Main output table</a></li>
				</ul>
			</div>
		</div>
	</div>

	<div class="pure-g"><div class="pure-u-1" id="divider_row"></div></div>
	<div class="pure-g" id="results_row">
		<div class="pure-u-2-24" id="margin"></div>
		<div class="pure-u-20-24" id="results_col">
			<div class="pure-g">
				<div class="pure-u-12-24">
					<h2>Yeast input phylogeny</h2>
					<img class="pure-img" id="logo_main" src="img/yeast_phylo.png">
					<h3>H1 Node: 20</h3>
					<h3>H2 Nodes: ZYGRO, TORDC, 9, 8, 22, 5, 7, 21, 20, 23, 24, HANAN</h3>
					<h3>505 trees after filtering</h3>
				</div>
				<div class="pure-u-12-24">
					<h2>Optimal Yeast MUL-tree</h2>
					<img class="pure-img" id="logo_main" src="img/yeast_mul_phylo.png">
					<h3>Optimal H1 node: 20</h3>
					<h3>Optimal H2 node: 22</h3>
				</div>
			</div>
			<div class="pure-g"><div class="pure-u-1" id="divider_row"></div></div>
			<div class="pure-g">
				<div class="pure-u-1">
					<h2>Score plot</h2>
					<center>
					{p1}
					</center>
				</div>
			</div>
			<div class="pure-g"><div class="pure-u-1" id="divider_row"></div></div>
			<div class="pure-g"><div class="pure-u-1" id="divider_row"></div></div>
		</div>
		<div class="pure-u-2-24" id="margin"></div>
	</div>

    {footer}
</body>
"""

######################
# Main block
######################
pagefile = "yeast_ohno.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Yeast Ohnologs"

head = RC.readResultsHead(title);
nav = RC.readYeastNav();
footer = RC.readFooter();

p1 = '<div id="bd28b00b-12ea-4318-9a4f-4c887ec3dff5" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("bd28b00b-12ea-4318-9a4f-4c887ec3dff5", [{"opacity": 0.6, "y": [25710, 26048, 26048, 26054, 26213, 26261, 26383, 26499, 26517, 26854, 26872, 26911, 28221], "type": "scatter", "mode": "markers", "x": ["<20>-<22>", "<20>-<21>", "<20>-<9>", "<20>-<23>", "<20>-HANAN", "<20>-<24>", "ST", "<20>-ZYGRO", "<20>-TORDC", "<20>-<8>", "<20>-<5>", "<20>-<7>", "<20>-<20>"]}], {"autosize": false, "title": "GRAMPA Results: yeast_ohno_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';

outfilename = "../../results/yeast/" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, p1=p1, footer=footer));