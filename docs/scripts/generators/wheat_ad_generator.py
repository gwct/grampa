############################################################
# For GRAMPA site, 12.19
# This generates the file "wheat_ad.html"
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

				<h1>Wheat Results  (AD sub-genomes)</h1>

				<div class="row">
					<a name="install"></a>
					<div class="col-1" id="jump_row">
						<div id="jump_container">
							<center>Files:
							<a class="jump_link" href="trees/wheat_A.tre" download>Species tree</a>
							<a class="jump_link" href="trees/wheat_AD_trees.txt" download>Gene trees</a>
							<a class="jump_link" href="output/wheat_AD_A_out.txt" download>Main output table</a>
							<a class="jump_link" href="wheat_all.html">All sub-genome results</a>
							<a class="jump_link" href="wheat_ab.html">AB sub-genome results</a>
							<a class="jump_link" href="wheat_bd.html">BD sub-genome results</a>
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
								<img class="grid-img" src="img/wheat_A_phylo.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>	
						<h3>H1 Node: All possible nodes</h3>
						<h3>H2 Nodes: All possible nodes</h3>
						<h3>7,958 trees after filtering</h3>									
					</div>

					<div class="col-12-24">
						<center><h2>Optimal Wheat MUL-tree</h2></center>
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/wheat_AD_mul.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>	
						<h3>Optimal H1 node: TAES (A sub-genome)</h3>
						<h3>Optimal H2 node: ATAU (D sub-genome)</h3>	
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
					<div class="col-1">
						<h2>Score plot with alternate species topology 1</h2>
						<center>
						{p2}
						</center>
					</div>
				</div>

				<div class="sep-div-2"></div>

				<center><h2>Known wheat phylogeny</h2></center>
				<div class="row img-row">							
					<div class="col-9-24 img-margin-left"></div>
					<div class="col-6-24 img-col">						
						<img class="grid-img" src="img/wheat_full_phylo.png">			
					</div>
					<div class="col-9-24 img-margin-right"></div>
				</div>	

				<div class="sep-div-2"></div>

				<div class="row">
					<div class="col-6-24 margin"></div>
					<div class="col-12-24">
						<h2>Species list</h2>
						<div class="table-container">
							<table class="stripe-table" id="spec_table">
								<thead><th>Species name</th><th>Species ID</th></thead>
								<tr><td>Setaria italica</td><td>SITA</td></tr>
								<tr><td>Zea mays</td><td>ZMAY</td></tr>
								<tr><td>Sorghum bicolor</td><td>SBIC</td></tr>
								<tr><td>Leersia perrieri</td><td>LPER</td></tr>
								<tr><td>Oryza sativa</td><td>OSAT</td></tr>
								<tr><td>Bracypodium distachyon</td><td>BDIS</td></tr>
								<tr><td>Hordeum vulgare</td><td>HVUL</td></tr>
								<tr><td>Triticum aestivum</td><td>TAES</td></tr>
								<tr><td>Triticum urartu</td><td>TURA</td></tr>
								<tr><td>Aegilops tauschii</td><td>ATAU</td></tr>
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
pagefile = "wheat_ad.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Wheat (AD)"

head = RC.readResultsHead(title);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

p1 = '<div id="94cd7898-db7f-4d73-b831-6da078e931ca" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("94cd7898-db7f-4d73-b831-6da078e931ca", [{"opacity": 0.6, "y": [75927, 82776, 82776, 90406, 93293, 94527, 94616, 94650, 95130, 95238, 95421, 95721, 96177, 96177, 96413, 96788, 96813, 97251, 97273, 97273, 98282, 98282, 98465, 98596, 98753, 98977, 99007, 99225, 99260, 99444, 99465, 99479, 99519, 99572, 99572, 99643, 99794, 99963, 100049, 100126, 100135, 100241, 100260, 100315, 100330, 100362, 100363, 100398, 100542, 100549, 100562, 100583, 100647, 100651, 100719, 100768, 100784, 100807, 100833, 100884, 100944, 101154, 101193, 101256, 101462, 101492, 101503, 101537, 101541, 101562, 101630, 101630, 101630, 101642, 101670, 101703, 101740, 101747, 101748, 101768, 101799, 101802, 101802, 101822, 101857, 101864, 101889, 101897, 101917, 101930, 101935, 101958, 101960, 101963, 102009, 102029, 102029, 102054, 102064, 102089, 102089, 102089, 102103, 102108, 102129, 102140, 102141, 102164, 102213, 102217, 102250, 102263, 102350, 102357, 102371, 102380, 102387, 102392, 102392, 102396, 102400, 102403, 102403, 102414, 102414, 102428, 102443, 102465, 102468, 102472, 102475, 102476, 102477, 102504, 102505, 102530, 102535, 102553, 102555, 102570, 102575, 102578, 102586, 102594, 102604, 102606, 102622, 102648, 102650, 102664, 102672, 102691, 102707, 102725, 102742, 102742, 102809, 102811, 102827, 102843, 102852, 102856, 102867, 102869, 102890, 102895, 102918, 102919, 102941, 102947, 102947, 103018, 103021, 103043, 103049, 103049, 103064, 103069, 103072, 103073, 103081, 103122, 103141, 103145, 103176, 103183, 103205, 103212, 103239, 103247, 103262, 103351, 103366, 103395, 103412, 103412, 103469, 103473, 103544, 103576, 103589, 103598, 103627, 103646, 103673, 103696, 103699, 103700, 103731, 103762, 103772, 103776, 103799, 103818, 103846, 103881, 103896, 103921, 103922, 103931, 104013, 104013, 104103, 104303, 104393, 104397, 104552, 104563, 104608, 104613, 104628, 104635, 104645, 104653, 104728, 104750, 104751, 104777, 104780, 104870, 104923, 104930, 104956, 104991, 105009, 105009, 105111, 105127, 105137, 105147, 105162, 105162, 105181, 105263, 105263, 105298, 105316, 105337, 105337, 105384, 105426, 105553, 105705, 105744, 105757, 105786, 105788, 105871, 105952, 106006, 106022, 106099, 106301, 106317, 106395, 106413, 106414, 106463, 106470, 106614, 107059, 108843, 109787, 109787, 111165, 111195, 111240, 111241, 111317, 111384, 111406, 111679, 111718, 111738, 112779, 112856, 112910, 112932, 113029, 113107, 113171, 113464, 113522, 113542], "type": "scatter", "mode": "markers", "x": ["TAES-ATAU", "<1>-ATAU", "<1>-<2>", "TAES-<2>", "ST", "<2>-<2>", "<2>-<9>", "<9>-<9>", "<3>-<9>", "<1>-<9>", "<1>-<3>", "<4>-<9>", "<6>-<8>", "<6>-<9>", "TAES-<9>", "<1>-HVUL", "TAES-HVUL", "TAES-<3>", "TAES-<1>", "TAES-TURA", "<2>-<3>", "<2>-HVUL", "ZMAY-ZMAY", "<1>-<4>", "<2>-<4>", "<7>-<9>", "TAES-BDIS", "<1>-BDIS", "TURA-<9>", "ZMAY-<9>", "<5>-<9>", "TURA-ATAU", "<1>-<5>", "<8>-<9>", "<8>-<6>", "ATAU-TURA", "TAES-<4>", "<2>-<6>", "SBIC-<9>", "TURA-HVUL", "TAES-<5>", "<1>-<6>", "<2>-<5>", "TAES-LPER", "OSAT-<9>", "<1>-LPER", "<1>-SITA", "<2>-BDIS", "<2>-<8>", "TAES-SITA", "<1>-<8>", "ATAU-<9>", "<1>-OSAT", "LPER-<9>", "HVUL-<9>", "TAES-OSAT", "SITA-<9>", "<2>-SITA", "BDIS-<9>", "TURA-<2>", "ATAU-TAES", "TURA-BDIS", "<3>-<6>", "TAES-<8>", "<3>-<8>", "<8>-HVUL", "TAES-<6>", "<7>-HVUL", "HVUL-HVUL", "ATAU-HVUL", "<3>-BDIS", "<3>-<4>", "<2>-LPER", "SBIC-HVUL", "<1>-<7>", "TURA-LPER", "SITA-SITA", "TAES-ZMAY", "ZMAY-HVUL", "BDIS-BDIS", "OSAT-HVUL", "ZMAY-SBIC", "ZMAY-<7>", "<5>-HVUL", "OSAT-OSAT", "<1>-ZMAY", "<2>-OSAT", "TAES-SBIC", "LPER-LPER", "TAES-<7>", "SBIC-SITA", "<3>-SITA", "<4>-<8>", "TURA-TURA", "<1>-SBIC", "SITA-HVUL", "LPER-HVUL", "<2>-<7>", "<8>-BDIS", "BDIS-HVUL", "OSAT-LPER", "OSAT-<5>", "ATAU-ATAU", "TURA-OSAT", "<3>-LPER", "TURA-<3>", "<3>-<5>", "ATAU-BDIS", "<7>-BDIS", "SBIC-SBIC", "TURA-SITA", "SBIC-BDIS", "<3>-OSAT", "<6>-SITA", "TURA-<5>", "ZMAY-BDIS", "<4>-SITA", "SBIC-ZMAY", "SBIC-<7>", "<8>-LPER", "<8>-ATAU", "OSAT-TURA", "<7>-ATAU", "<7>-SITA", "<7>-<8>", "<2>-ZMAY", "<8>-<5>", "<7>-TURA", "<4>-LPER", "OSAT-BDIS", "<5>-TURA", "<8>-TURA", "HVUL-BDIS", "SBIC-ATAU", "<5>-BDIS", "OSAT-ATAU", "<7>-LPER", "ZMAY-SITA", "<2>-SBIC", "<5>-ATAU", "ATAU-LPER", "<3>-ZMAY", "SBIC-LPER", "ZMAY-ATAU", "SBIC-TURA", "LPER-BDIS", "HVUL-TURA", "<7>-<7>", "ZMAY-TURA", "SITA-BDIS", "<3>-SBIC", "ZMAY-LPER", "HVUL-ATAU", "SBIC-<8>", "LPER-OSAT", "LPER-<5>", "SITA-LPER", "LPER-ATAU", "<3>-<3>", "BDIS-LPER", "SITA-ATAU", "HVUL-LPER", "<4>-OSAT", "BDIS-ATAU", "LPER-TURA", "SITA-TURA", "<8>-OSAT", "ZMAY-<8>", "<4>-ZMAY", "<5>-SITA", "BDIS-TURA", "<4>-SBIC", "<7>-OSAT", "<6>-ZMAY", "<4>-<6>", "<4>-<5>", "<3>-<7>", "<6>-SBIC", "ATAU-SITA", "OSAT-SITA", "TURA-ZMAY", "SBIC-OSAT", "ATAU-OSAT", "HVUL-SITA", "<5>-<8>", "LPER-SITA", "TURA-SBIC", "ZMAY-OSAT", "<5>-<5>", "BDIS-SITA", "TURA-<4>", "<7>-<5>", "SITA-OSAT", "<4>-<7>", "HVUL-OSAT", "BDIS-OSAT", "<6>-<7>", "SBIC-<5>", "ATAU-<5>", "BDIS-<5>", "<5>-ZMAY", "TURA-<8>", "ZMAY-<5>", "OSAT-ZMAY", "SITA-ZMAY", "SITA-SBIC", "<5>-SBIC", "ATAU-ZMAY", "HVUL-ZMAY", "TURA-<7>", "OSAT-SBIC", "LPER-ZMAY", "BDIS-ZMAY", "ATAU-SBIC", "HVUL-SBIC", "HVUL-<5>", "LPER-SBIC", "SITA-<5>", "BDIS-SBIC", "<8>-<4>", "SITA-<7>", "SITA-<8>", "<8>-<8>", "ATAU-<3>", "<5>-<7>", "OSAT-<8>", "OSAT-<4>", "OSAT-<7>", "ATAU-<7>", "LPER-<8>", "<8>-<3>", "<7>-<3>", "ATAU-<8>", "HVUL-<7>", "LPER-<7>", "<7>-<4>", "HVUL-<8>", "BDIS-<8>", "BDIS-<7>", "TURA-<6>", "OSAT-<3>", "<7>-<6>", "ATAU-<4>", "SBIC-<3>", "HVUL-<2>", "HVUL-<3>", "ZMAY-<3>", "<4>-<4>", "SBIC-<4>", "<5>-<3>", "BDIS-<4>", "BDIS-<3>", "LPER-<4>", "<5>-<6>", "<5>-<4>", "ZMAY-<4>", "OSAT-<6>", "ATAU-<1>", "ATAU-<2>", "HVUL-<4>", "LPER-<3>", "SITA-<3>", "<7>-<2>", "SITA-<4>", "<8>-<2>", "ZMAY-<6>", "SBIC-<6>", "OSAT-<2>", "LPER-<6>", "SBIC-<2>", "<5>-<2>", "ZMAY-<2>", "BDIS-<6>", "ATAU-<6>", "LPER-<2>", "BDIS-<2>", "<6>-<6>", "SITA-<6>", "SITA-<2>", "HVUL-<6>", "TAES-TAES", "<1>-<1>", "TURA-TAES", "TURA-<1>", "<7>-TAES", "HVUL-TAES", "OSAT-TAES", "<8>-TAES", "<5>-TAES", "SBIC-TAES", "ZMAY-TAES", "LPER-TAES", "SITA-TAES", "BDIS-TAES", "HVUL-<1>", "<7>-<1>", "OSAT-<1>", "<8>-<1>", "<5>-<1>", "SBIC-<1>", "ZMAY-<1>", "LPER-<1>", "SITA-<1>", "BDIS-<1>"]}], {"autosize": false, "title": "GRAMPA Results: wheat_ad_A_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';
p2 = '<div id="6adfb019-e879-438f-9c76-1c331069c044" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("6adfb019-e879-438f-9c76-1c331069c044", [{"opacity": 0.6, "y": [75955, 82684, 82684, 90089, 92946, 94181, 94242, 94296, 94765, 95354, 95379, 95656, 95845, 95845, 96049, 96397, 96435, 96671, 96671, 96866, 97923, 97923, 98114, 98405, 98643, 98675, 98887, 99037, 99102, 99138, 99218, 99222, 99227, 99227, 99407, 99437, 99604, 99721, 99766, 99769, 99912, 99971, 99979, 99995, 100070, 100079, 100190, 100198, 100305, 100359, 100381, 100423, 100436, 100444, 100465, 100490, 100566, 100698, 100846, 100880, 100884, 100952, 101020, 101131, 101132, 101155, 101184, 101185, 101191, 101195, 101297, 101299, 101299, 101306, 101391, 101400, 101411, 101438, 101454, 101460, 101460, 101479, 101524, 101535, 101549, 101562, 101573, 101602, 101613, 101630, 101632, 101684, 101686, 101689, 101708, 101737, 101745, 101746, 101747, 101751, 101757, 101757, 101793, 101800, 101866, 101885, 101890, 101927, 101938, 101945, 102009, 102039, 102048, 102058, 102058, 102059, 102063, 102064, 102071, 102078, 102078, 102094, 102099, 102099, 102132, 102133, 102135, 102145, 102154, 102176, 102179, 102205, 102208, 102217, 102221, 102228, 102230, 102251, 102254, 102259, 102267, 102268, 102269, 102276, 102305, 102327, 102336, 102340, 102344, 102355, 102363, 102370, 102393, 102393, 102407, 102454, 102482, 102490, 102498, 102501, 102504, 102513, 102524, 102527, 102528, 102530, 102533, 102572, 102580, 102582, 102618, 102620, 102682, 102684, 102694, 102701, 102701, 102717, 102734, 102738, 102739, 102794, 102814, 102829, 102851, 102858, 102864, 102896, 102914, 103020, 103032, 103072, 103074, 103075, 103095, 103145, 103154, 103207, 103208, 103237, 103254, 103274, 103305, 103316, 103341, 103355, 103358, 103397, 103423, 103444, 103467, 103502, 103545, 103553, 103579, 103587, 103623, 103668, 103668, 103746, 103759, 103789, 103792, 103849, 104054, 104055, 104055, 104063, 104086, 104216, 104231, 104276, 104313, 104315, 104316, 104389, 104413, 104437, 104438, 104440, 104593, 104595, 104676, 104688, 104688, 104778, 104798, 104817, 104833, 104833, 104834, 104851, 104880, 104880, 104932, 104932, 104961, 104985, 105060, 105094, 105221, 105312, 105388, 105417, 105428, 105445, 105472, 105479, 105551, 105610, 105697, 105703, 105771, 105955, 106065, 106072, 106087, 106124, 106146, 106279, 106725, 108196, 109486, 109486, 110839, 110891, 110939, 111017, 111073, 111076, 111081, 111344, 111393, 111394, 112218, 112332, 112395, 112513, 112589, 112602, 112619, 112905, 112952, 112965], "type": "scatter", "mode": "markers", "x": ["TAES-TURA", "<1>-<2>", "<1>-TURA", "TAES-<2>", "ST", "<2>-<2>", "<2>-<9>", "<9>-<9>", "<3>-<9>", "<1>-<3>", "<4>-<9>", "<1>-<9>", "<6>-<8>", "<6>-<9>", "TAES-<9>", "<1>-HVUL", "TAES-HVUL", "TAES-ATAU", "TAES-<1>", "TAES-<3>", "<2>-<3>", "<2>-HVUL", "ZMAY-ZMAY", "<2>-<4>", "<7>-<9>", "TAES-BDIS", "<1>-<4>", "<1>-BDIS", "ZMAY-<9>", "<5>-<9>", "ATAU-TURA", "TURA-ATAU", "<8>-<9>", "<8>-<6>", "TURA-<9>", "TAES-<4>", "<2>-<6>", "SBIC-<9>", "<1>-<5>", "TAES-<5>", "<2>-<5>", "ATAU-<9>", "TAES-LPER", "OSAT-<9>", "ATAU-HVUL", "<2>-BDIS", "TAES-SITA", "<2>-<8>", "LPER-<9>", "<1>-SITA", "HVUL-<9>", "TAES-OSAT", "SITA-<9>", "<1>-LPER", "<2>-SITA", "BDIS-<9>", "<1>-<6>", "<1>-<8>", "<3>-<6>", "TAES-<8>", "<1>-OSAT", "TURA-TAES", "ATAU-<2>", "<3>-<8>", "TAES-<6>", "<8>-HVUL", "ATAU-BDIS", "TURA-HVUL", "HVUL-HVUL", "<7>-HVUL", "<2>-LPER", "<3>-BDIS", "<3>-<4>", "SBIC-HVUL", "TAES-ZMAY", "ZMAY-HVUL", "SITA-SITA", "BDIS-BDIS", "OSAT-HVUL", "ZMAY-SBIC", "ZMAY-<7>", "<5>-HVUL", "OSAT-OSAT", "TAES-SBIC", "<2>-OSAT", "TAES-<7>", "LPER-LPER", "SBIC-SITA", "TURA-TURA", "<3>-SITA", "<4>-<8>", "LPER-HVUL", "SITA-HVUL", "<1>-<7>", "<2>-<7>", "ATAU-ATAU", "<8>-BDIS", "<1>-ZMAY", "BDIS-HVUL", "TURA-BDIS", "OSAT-LPER", "OSAT-<5>", "<3>-LPER", "<3>-<5>", "SBIC-SBIC", "<1>-SBIC", "<7>-BDIS", "TURA-LPER", "ATAU-LPER", "SBIC-BDIS", "<3>-OSAT", "<6>-SITA", "ZMAY-BDIS", "SBIC-ZMAY", "SBIC-<7>", "<8>-TURA", "<7>-TURA", "<4>-SITA", "<8>-LPER", "<7>-SITA", "<7>-<8>", "<2>-ZMAY", "OSAT-TURA", "<7>-ATAU", "<8>-ATAU", "<8>-<5>", "<4>-LPER", "OSAT-BDIS", "HVUL-BDIS", "<5>-TURA", "<5>-BDIS", "SBIC-TURA", "<7>-LPER", "<2>-SBIC", "ZMAY-SITA", "SBIC-ATAU", "OSAT-ATAU", "<3>-ZMAY", "<5>-ATAU", "ZMAY-TURA", "SBIC-LPER", "HVUL-ATAU", "ZMAY-ATAU", "LPER-BDIS", "<7>-<7>", "TURA-OSAT", "SITA-BDIS", "<3>-SBIC", "TURA-SITA", "ZMAY-LPER", "HVUL-TURA", "ATAU-<3>", "LPER-OSAT", "LPER-<5>", "SBIC-<8>", "LPER-ATAU", "SITA-LPER", "<3>-<3>", "ATAU-OSAT", "BDIS-ATAU", "SITA-ATAU", "BDIS-LPER", "ATAU-SITA", "<4>-OSAT", "HVUL-LPER", "LPER-TURA", "SITA-TURA", "ZMAY-<8>", "BDIS-TURA", "<8>-OSAT", "<5>-SITA", "<4>-ZMAY", "<7>-OSAT", "TURA-<5>", "<4>-SBIC", "<4>-<6>", "<4>-<5>", "<6>-ZMAY", "<3>-<7>", "OSAT-SITA", "<6>-SBIC", "SBIC-OSAT", "HVUL-SITA", "ATAU-<5>", "LPER-SITA", "<5>-<8>", "ZMAY-OSAT", "<5>-<5>", "BDIS-SITA", "<7>-<5>", "SITA-OSAT", "<4>-<7>", "HVUL-OSAT", "BDIS-OSAT", "TURA-ZMAY", "<6>-<7>", "SBIC-<5>", "ATAU-ZMAY", "TURA-SBIC", "BDIS-<5>", "<5>-ZMAY", "ZMAY-<5>", "OSAT-ZMAY", "ATAU-SBIC", "SITA-ZMAY", "SITA-SBIC", "<5>-SBIC", "HVUL-ZMAY", "OSAT-SBIC", "LPER-ZMAY", "BDIS-ZMAY", "HVUL-SBIC", "HVUL-<5>", "LPER-SBIC", "BDIS-SBIC", "SITA-<5>", "<8>-<4>", "SITA-<7>", "SITA-<8>", "TURA-<8>", "<8>-<8>", "ATAU-<4>", "TURA-<3>", "TURA-<7>", "<5>-<7>", "OSAT-<8>", "ATAU-<7>", "ATAU-<8>", "TURA-<4>", "OSAT-<7>", "OSAT-<4>", "LPER-<8>", "HVUL-<7>", "<8>-<3>", "<7>-<3>", "LPER-<7>", "HVUL-<8>", "BDIS-<8>", "<7>-<4>", "BDIS-<7>", "<7>-<6>", "OSAT-<3>", "SBIC-<3>", "HVUL-<2>", "HVUL-<3>", "ZMAY-<3>", "<4>-<4>", "<5>-<3>", "BDIS-<4>", "BDIS-<3>", "SBIC-<4>", "LPER-<4>", "TURA-<2>", "TURA-<1>", "<5>-<6>", "<5>-<4>", "ZMAY-<4>", "OSAT-<6>", "HVUL-<4>", "LPER-<3>", "SITA-<3>", "TURA-<6>", "<7>-<2>", "SITA-<4>", "ZMAY-<6>", "<8>-<2>", "SBIC-<6>", "ATAU-<6>", "OSAT-<2>", "LPER-<6>", "SBIC-<2>", "<5>-<2>", "ZMAY-<2>", "BDIS-<6>", "<6>-<6>", "LPER-<2>", "BDIS-<2>", "SITA-<6>", "SITA-<2>", "HVUL-<6>", "TAES-TAES", "<1>-<1>", "ATAU-<1>", "ATAU-TAES", "<7>-TAES", "HVUL-TAES", "<8>-TAES", "OSAT-TAES", "SBIC-TAES", "<5>-TAES", "ZMAY-TAES", "LPER-TAES", "SITA-TAES", "BDIS-TAES", "HVUL-<1>", "<7>-<1>", "<8>-<1>", "OSAT-<1>", "SBIC-<1>", "<5>-<1>", "ZMAY-<1>", "LPER-<1>", "BDIS-<1>", "SITA-<1>"]}], {"autosize": false, "title": "GRAMPA Results: wheat_ad_D_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';

outfilename = "../../results/wheat/" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, p1=p1, p2=p2, footer=footer));