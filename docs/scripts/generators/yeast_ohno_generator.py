############################################################
# For GRAMPA site, 12.19
# This generates the file "yeast_ohno.html"
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
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_col">
			<div id="main_content">

				<h1>Yeast Results</h1>

				<div class="row">
					<a name="install"></a>
					<div class="col-1" id="jump_row">
						<div id="jump_container">
							<center>Files:
							<a class="jump_link" href="trees/yeast_spec.tre" download>Species tree</a>
							<a class="jump_link" href="trees/yeast_ohno_trees_final.txt" download>Gene trees</a>
							<a class="jump_link" href="output/yeast_ohno_out.txt" download>Main output table</a>
							<a class="jump_link" href="yeast.html">All homolog results</a>
							</center>
						</div>
					</div>
				</div>					

				<div class="row">
					<div class="col-12-24">
						<center><h2>Input phylogeny</h2></center>
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/yeast_phylo.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>	
						<h3>H1 Node: 20</h3>
						<h3>H2 Nodes: ZYGRO, TORDC, 9, 8, 22, 5, 7, 21, 20, 23, 24, HANAN</h3>									
					</div>

					<div class="col-12-24">
						<center><h2>Optimal Yeast MUL-tree</h2></center>
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/yeast_mul_phylo.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>	
						<h3>Optimal H1 node: 20</h3>
						<h3>Optimal H2 node: 22</h3>					
					</div>
				</div>

				<div class="sep-div-2"></div>
				
				<div class="row">
					<div class="col-1">
						<h2>Score plot</h2>
						<center>
						{p1}
						</center>
					</div>
				</div>

				<div class="sep-div-2"></div>

				<div class="row">
					<div class="col-6-24 margin"></div>
					<div class="col-12-24">
						<h2>Species list</h2>
						<div class="table-container">
							<table class="stripe-table" id="spec_table">
								<thead><th>Species name</th><th>AKA</th><th>Species ID</th></thead>
								<tr><td>Schizosaccharomyces pombe</td><td></td><td>SCHPO</td></tr>
								<tr><td>Yarrowia lipolytica</td><td></td><td>YARLI</td></tr>
								<tr><td>Nadsonia fulvescens</td><td></td><td>45786A</td></tr>
								<tr><td>Dekkera bruxellensis</td><td></td><td>DEKBR</td></tr>
								<tr><td>Debaryomyces hansenii</td><td></td><td>DEBHA</td></tr>
								<tr><td>Scheffersomyces stipitis</td><td>Pichia stipitis</td><td>PICST</td></tr>
								<tr><td>Candida albicans</td><td></td><td>CANAL</td></tr>
								<tr><td>Wickerhamomyces anomalus</td><td></td><td>HANAN</td></tr>
								<tr><td>Kluyveromyces lactis</td><td></td><td>KLULA</td></tr>
								<tr><td>Ashbya gossypii</td><td></td><td>ASHGO</td></tr>
								<tr><td>Lachancea kluyveri</td><td>Saccharomyces kluyveri</td><td>SACKL</td></tr>
								<tr><td>Lachancea waltii</td><td>Kluyveromyces waltii</td><td>KLUWA</td></tr>
								<tr><td>Lachancea thermotolerans</td><td></td><td>LACTH</td></tr>
								<tr><td>Torulaspora delbrueckii</td><td></td><td>TORDC</td></tr>
								<tr><td>Zygosaccharomyces rouxii</td><td></td><td>ZYGRO</td></tr>
								<tr><td>Tetrapisispora blattae</td><td></td><td>1071379A</td></tr>
								<tr><td>Tetrapisispora phaffii</td><td></td><td>TETPH</td></tr>
								<tr><td>Vanderwaltozyma polyspora</td><td></td><td>VANPO</td></tr>
								<tr><td>Kazachstania naganishii</td><td></td><td>588726A</td></tr>
								<tr><td>Kazachstania africana</td><td></td><td>KAZAF</td></tr>
								<tr><td>Naumovozyma dairenensis</td><td></td><td>NAUDC</td></tr>
								<tr><td>Naumovozyma castellii</td><td>Saccharomyces castellii</td><td>SACCA</td></tr>
								<tr><td>Saccharomyces bayanus</td><td></td><td>SACBA</td></tr>
								<tr><td>Saccharomyces cerevisiae</td><td></td><td>YEAST</td></tr>
								<tr><td>Candida glabrata</td><td></td><td>CANGA</td></tr>
								<tr><td>Nakaseomyces bacillisporus</td><td></td><td>51600A</td></tr>
								<tr><td>Candida castellii</td><td></td><td>51914A</td></tr>
							</table>
						</div>
					</div>
					<div class="col-6-24 margin"></div>
				</div>
			</div>
			<div class="col-2-24" id="margin"></div>
		</div>
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
nav = RC.readNav(pagefile, "../../", "", "../wheat/");
footer = RC.readFooter();

p1 = '<div id="bd28b00b-12ea-4318-9a4f-4c887ec3dff5" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("bd28b00b-12ea-4318-9a4f-4c887ec3dff5", [{"opacity": 0.6, "y": [25710, 26048, 26048, 26054, 26213, 26261, 26383, 26499, 26517, 26854, 26872, 26911, 28221], "type": "scatter", "mode": "markers", "x": ["<20>-<22>", "<20>-<21>", "<20>-<9>", "<20>-<23>", "<20>-HANAN", "<20>-<24>", "ST", "<20>-ZYGRO", "<20>-TORDC", "<20>-<8>", "<20>-<5>", "<20>-<7>", "<20>-<20>"]}], {"autosize": false, "title": "GRAMPA Results: yeast_ohno_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';

outfilename = "../../results/yeast/" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, p1=p1, footer=footer));