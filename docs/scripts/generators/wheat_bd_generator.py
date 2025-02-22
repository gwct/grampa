############################################################
# For GRAMPA site, 12.19
# This generates the file "wheat_bd.html"
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

				<h1>Wheat Results  (BD sub-genomes)</h1>

				<div class="row">
					<a name="install"></a>
					<div class="col-1" id="jump_row">
						<div id="jump_container">
							<center>Files:
							<a class="jump_link" href="trees/wheat_D.tre" download>Species tree</a>
							<a class="jump_link" href="trees/wheat_BD_trees.txt" download>Gene trees</a>
							<a class="jump_link" href="output/wheat_BD_D_out.txt" download>Main output table</a>
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
						<h3>7,986 trees after filtering</h3>									
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
						<h3>Optimal H1 node: ATAU (D sub-genome)</h3>
						<h3>Optimal H2 node: 2 (B sub-genome)</h3>
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
pagefile = "wheat_bd.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Wheat (BD)"

head = RC.readResultsHead(title);
nav = RC.readNav(pagefile, "../../", "../yeast/", "");
footer = RC.readFooter();

p1 = '<div id="52046d92-80e0-447e-8f2c-d069bda5f982" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("52046d92-80e0-447e-8f2c-d069bda5f982", [{"opacity": 0.6, "y": [82796, 85682, 85682, 87987, 91187, 91368, 91368, 92273, 92508, 92862, 93532, 93663, 93783, 93826, 93944, 93969, 94049, 94049, 94836, 95419, 96266, 96266, 96266, 96318, 96500, 96501, 96546, 96934, 96961, 97140, 97149, 97378, 97384, 97509, 97509, 97673, 97673, 97783, 97925, 97995, 98059, 98077, 98203, 98239, 98280, 98297, 98316, 98487, 98525, 98575, 98628, 98635, 98642, 98672, 98686, 98690, 98757, 98792, 98878, 98995, 99050, 99086, 99288, 99374, 99412, 99423, 99434, 99475, 99495, 99500, 99500, 99516, 99582, 99605, 99639, 99639, 99685, 99686, 99690, 99699, 99711, 99726, 99726, 99762, 99763, 99804, 99818, 99840, 99851, 99856, 99861, 99891, 99907, 99907, 99981, 99987, 99987, 100007, 100014, 100030, 100038, 100052, 100069, 100124, 100130, 100180, 100190, 100239, 100275, 100277, 100280, 100284, 100287, 100308, 100309, 100309, 100315, 100322, 100322, 100337, 100338, 100367, 100367, 100378, 100389, 100402, 100410, 100418, 100423, 100431, 100450, 100471, 100477, 100479, 100492, 100515, 100516, 100520, 100534, 100548, 100551, 100553, 100556, 100563, 100565, 100569, 100586, 100614, 100614, 100617, 100645, 100648, 100729, 100742, 100745, 100746, 100754, 100754, 100776, 100781, 100786, 100790, 100804, 100820, 100831, 100839, 100840, 100891, 100919, 100929, 100930, 100936, 100936, 100964, 100976, 100986, 100998, 101004, 101021, 101085, 101101, 101101, 101117, 101129, 101185, 101207, 101266, 101286, 101287, 101296, 101306, 101307, 101313, 101376, 101401, 101402, 101443, 101503, 101512, 101513, 101534, 101554, 101572, 101584, 101612, 101624, 101656, 101675, 101681, 101713, 101766, 101805, 101827, 101839, 101841, 101911, 101911, 101935, 101986, 102050, 102155, 102158, 102287, 102320, 102333, 102339, 102369, 102468, 102470, 102506, 102528, 102583, 102622, 102627, 102634, 102647, 102690, 102701, 102711, 102714, 102885, 102960, 103006, 103068, 103114, 103115, 103123, 103150, 103150, 103180, 103180, 103181, 103225, 103225, 103252, 103273, 103427, 103453, 103453, 103467, 103566, 103670, 103715, 103736, 103773, 103795, 103858, 103882, 103936, 104025, 104124, 104163, 104217, 104289, 104336, 104411, 104454, 104498, 104529, 104591, 104613, 104882, 106655, 106655, 108959, 109112, 109169, 109186, 109239, 109268, 109355, 109386, 109422, 109423, 109436, 109443, 109446, 109477, 109693, 109707, 109742, 109752, 109756, 109761], "type": "scatter", "mode": "markers", "x": ["TAES-<2>", "<1>-<2>", "<1>-TURA", "TAES-TURA", "ST", "TAES-ATAU", "TAES-<1>", "<2>-<9>", "<9>-<9>", "<3>-<9>", "<4>-<9>", "<1>-<9>", "TAES-<3>", "TAES-HVUL", "<1>-<3>", "TAES-<9>", "<6>-<8>", "<6>-<9>", "<1>-HVUL", "<2>-<2>", "<2>-<3>", "ATAU-TURA", "<2>-HVUL", "ZMAY-ZMAY", "TAES-BDIS", "<2>-<4>", "TURA-ATAU", "<7>-<9>", "<1>-<4>", "<1>-BDIS", "TAES-<4>", "ZMAY-<9>", "<5>-<9>", "<8>-<9>", "<8>-<6>", "TURA-<9>", "<2>-<6>", "TAES-<5>", "<1>-<5>", "SBIC-<9>", "<2>-<5>", "TAES-LPER", "OSAT-<9>", "TAES-SITA", "<2>-BDIS", "<2>-<8>", "ATAU-<9>", "<1>-SITA", "TAES-OSAT", "LPER-<9>", "<2>-SITA", "<1>-LPER", "HVUL-<9>", "<1>-<6>", "ATAU-HVUL", "SITA-<9>", "BDIS-<9>", "<1>-<8>", "TAES-<8>", "<3>-<6>", "TAES-<6>", "<1>-OSAT", "<3>-<8>", "<8>-HVUL", "<7>-HVUL", "HVUL-HVUL", "TURA-HVUL", "TAES-ZMAY", "<2>-LPER", "<3>-BDIS", "<3>-<4>", "SBIC-HVUL", "ATAU-BDIS", "TAES-<7>", "TAES-SBIC", "ZMAY-HVUL", "TURA-TURA", "SITA-SITA", "OSAT-HVUL", "BDIS-BDIS", "<5>-HVUL", "ZMAY-SBIC", "ZMAY-<7>", "OSAT-OSAT", "<2>-OSAT", "<1>-<7>", "LPER-LPER", "<3>-SITA", "<2>-<7>", "SBIC-SITA", "<4>-<8>", "<1>-ZMAY", "SITA-HVUL", "LPER-HVUL", "BDIS-HVUL", "OSAT-LPER", "OSAT-<5>", "<8>-BDIS", "<3>-LPER", "TURA-BDIS", "<3>-<5>", "<1>-SBIC", "ATAU-ATAU", "<7>-BDIS", "SBIC-SBIC", "TURA-LPER", "SBIC-BDIS", "<3>-OSAT", "<2>-ZMAY", "<8>-TURA", "ATAU-LPER", "<7>-TURA", "ZMAY-BDIS", "<4>-SITA", "SBIC-ZMAY", "SBIC-<7>", "<6>-SITA", "<7>-SITA", "<7>-<8>", "<8>-LPER", "OSAT-TURA", "<7>-ATAU", "<4>-LPER", "<8>-ATAU", "OSAT-BDIS", "<5>-TURA", "<8>-<5>", "<2>-SBIC", "SBIC-TURA", "<5>-BDIS", "<3>-ZMAY", "<7>-LPER", "HVUL-BDIS", "ZMAY-SITA", "ZMAY-TURA", "SBIC-LPER", "LPER-BDIS", "SBIC-ATAU", "HVUL-ATAU", "ZMAY-ATAU", "HVUL-TURA", "<3>-SBIC", "OSAT-ATAU", "TURA-OSAT", "<5>-ATAU", "<7>-<7>", "SITA-BDIS", "LPER-OSAT", "LPER-<5>", "ZMAY-LPER", "SBIC-<8>", "TURA-SITA", "SITA-LPER", "LPER-ATAU", "LPER-TURA", "<4>-OSAT", "SITA-TURA", "BDIS-LPER", "HVUL-LPER", "SITA-ATAU", "BDIS-ATAU", "BDIS-TURA", "ATAU-SITA", "ATAU-OSAT", "ZMAY-<8>", "<8>-OSAT", "<4>-ZMAY", "<5>-SITA", "<3>-<7>", "<4>-SBIC", "<7>-OSAT", "<4>-<6>", "<4>-<5>", "TURA-<5>", "<3>-<3>", "<6>-ZMAY", "OSAT-SITA", "<6>-SBIC", "SBIC-OSAT", "HVUL-SITA", "LPER-SITA", "<5>-<8>", "ZMAY-OSAT", "<5>-<5>", "BDIS-SITA", "ATAU-<5>", "SITA-OSAT", "ATAU-<2>", "<4>-<7>", "<7>-<5>", "HVUL-OSAT", "BDIS-OSAT", "ATAU-<3>", "TURA-ZMAY", "SBIC-<5>", "<6>-<7>", "ATAU-ZMAY", "TURA-SBIC", "BDIS-<5>", "<5>-ZMAY", "ZMAY-<5>", "OSAT-ZMAY", "ATAU-SBIC", "SITA-ZMAY", "SITA-SBIC", "<5>-SBIC", "HVUL-ZMAY", "OSAT-SBIC", "LPER-ZMAY", "BDIS-ZMAY", "HVUL-SBIC", "LPER-SBIC", "HVUL-<5>", "BDIS-SBIC", "SITA-<5>", "SITA-<7>", "SITA-<8>", "<8>-<4>", "<8>-<8>", "TURA-<8>", "TURA-<7>", "TURA-<3>", "OSAT-<8>", "<5>-<7>", "ATAU-<7>", "ATAU-<4>", "ATAU-<8>", "OSAT-<7>", "TURA-<4>", "LPER-<8>", "OSAT-<4>", "HVUL-<7>", "<7>-<3>", "<1>-<1>", "LPER-<7>", "<8>-<3>", "HVUL-<8>", "BDIS-<7>", "<7>-<4>", "BDIS-<8>", "<7>-<6>", "OSAT-<3>", "SBIC-<3>", "<4>-<4>", "SBIC-<4>", "ZMAY-<3>", "LPER-<4>", "HVUL-<2>", "HVUL-<3>", "BDIS-<4>", "BDIS-<3>", "<5>-<3>", "<5>-<6>", "<5>-<4>", "ZMAY-<4>", "OSAT-<6>", "LPER-<3>", "TURA-<2>", "TURA-<1>", "HVUL-<4>", "SITA-<3>", "TURA-<6>", "SITA-<4>", "ZMAY-<6>", "SBIC-<6>", "<7>-<2>", "<8>-<2>", "LPER-<6>", "ATAU-<6>", "OSAT-<2>", "SBIC-<2>", "<5>-<2>", "ZMAY-<2>", "BDIS-<6>", "<6>-<6>", "SITA-<6>", "TAES-TAES", "LPER-<2>", "BDIS-<2>", "SITA-<2>", "HVUL-<6>", "TURA-TAES", "ATAU-<1>", "ATAU-TAES", "HVUL-<1>", "<7>-<1>", "<8>-<1>", "<7>-TAES", "HVUL-TAES", "<8>-TAES", "OSAT-<1>", "SBIC-<1>", "ZMAY-<1>", "<5>-<1>", "OSAT-TAES", "SBIC-TAES", "ZMAY-TAES", "<5>-TAES", "LPER-<1>", "LPER-TAES", "BDIS-<1>", "SITA-TAES", "BDIS-TAES", "SITA-<1>"]}], {"autosize": false, "title": "GRAMPA Results: wheat_bd_D_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';
p2 = '<div id="a644023f-5e5c-4d74-9452-b50f7c186baf" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("a644023f-5e5c-4d74-9452-b50f7c186baf", [{"opacity": 0.6, "y": [83878, 96036, 96679, 96679, 97583, 97837, 100263, 100819, 100819, 101306, 101531, 101906, 102582, 102863, 103108, 103108, 104141, 105300, 105300, 105424, 105494, 105625, 105668, 106013, 106160, 106277, 106453, 106482, 106536, 106608, 106608, 106857, 106867, 107099, 107243, 107310, 107345, 107407, 107419, 107478, 107688, 107728, 107741, 107746, 107800, 107804, 107861, 107891, 108140, 108199, 108350, 108409, 108438, 108483, 108487, 108555, 108593, 108620, 108646, 108665, 108688, 108695, 108695, 108704, 108764, 108792, 108797, 108822, 108883, 108884, 108893, 108894, 108894, 108907, 108920, 108951, 108953, 108966, 109015, 109027, 109038, 109056, 109059, 109070, 109100, 109104, 109119, 109175, 109180, 109188, 109192, 109192, 109192, 109209, 109214, 109218, 109218, 109237, 109278, 109332, 109334, 109354, 109354, 109405, 109433, 109440, 109463, 109477, 109481, 109494, 109506, 109509, 109530, 109530, 109533, 109541, 109541, 109550, 109563, 109579, 109599, 109607, 109609, 109617, 109619, 109637, 109650, 109655, 109656, 109662, 109664, 109683, 109700, 109701, 109714, 109717, 109738, 109740, 109742, 109745, 109758, 109772, 109794, 109805, 109806, 109807, 109813, 109824, 109840, 109840, 109840, 109867, 109878, 109931, 109952, 109955, 109956, 109962, 109963, 109970, 109980, 110003, 110006, 110009, 110018, 110031, 110051, 110063, 110093, 110102, 110122, 110131, 110140, 110147, 110155, 110155, 110204, 110210, 110216, 110227, 110236, 110278, 110294, 110317, 110317, 110329, 110329, 110353, 110379, 110386, 110400, 110401, 110484, 110502, 110506, 110512, 110522, 110526, 110608, 110618, 110622, 110672, 110682, 110733, 110747, 110749, 110775, 110792, 110794, 110823, 110831, 110834, 110885, 110890, 110894, 110921, 110953, 110961, 110979, 111020, 111052, 111054, 111064, 111153, 111153, 111167, 111236, 111536, 111552, 111556, 111595, 111596, 111708, 111753, 111765, 111818, 111862, 111871, 111895, 111932, 111939, 111947, 111953, 112086, 112131, 112211, 112260, 112300, 112359, 112366, 112380, 112406, 112406, 112429, 112429, 112442, 112484, 112484, 112491, 112530, 112694, 112717, 112827, 112970, 112971, 113035, 113036, 113103, 113151, 113188, 113280, 113376, 113424, 113464, 113552, 113598, 113604, 113670, 113678, 113763, 113789, 113852, 113879, 116696, 117138, 117139, 117175, 117236, 117247, 117290, 117368, 117574, 117605, 117619, 118483, 118531, 118604, 118707, 118712, 118755, 118766, 119042, 119069, 119099], "type": "scatter", "mode": "markers", "x": ["TAES-ATAU", "TAES-TURA", "TAES-<1>", "TAES-<2>", "ATAU-TAES", "TURA-<2>", "ST", "<1>-<2>", "<1>-TAES", "<2>-<9>", "<9>-<9>", "<3>-<9>", "<4>-<9>", "<2>-<2>", "<6>-<8>", "<6>-<9>", "TAES-<9>", "<2>-<3>", "<2>-HVUL", "<1>-<9>", "ZMAY-ZMAY", "TAES-HVUL", "<2>-<4>", "<7>-<9>", "TURA-HVUL", "TURA-<9>", "ZMAY-<9>", "<5>-<9>", "<1>-HVUL", "<8>-<9>", "<8>-<6>", "TAES-<3>", "<2>-<6>", "SBIC-<9>", "<2>-<5>", "OSAT-<9>", "<1>-<3>", "<2>-BDIS", "ATAU-<9>", "<2>-<8>", "LPER-<9>", "<1>-BDIS", "TAES-BDIS", "HVUL-<9>", "<2>-SITA", "SITA-<9>", "BDIS-<9>", "ATAU-HVUL", "TURA-BDIS", "<3>-<6>", "TURA-<3>", "TAES-SITA", "TAES-<5>", "TAES-LPER", "<3>-<8>", "<8>-HVUL", "<7>-HVUL", "HVUL-HVUL", "TAES-<4>", "<1>-LPER", "<2>-LPER", "<3>-BDIS", "<3>-<4>", "SBIC-HVUL", "<1>-<5>", "ATAU-BDIS", "TURA-LPER", "ZMAY-HVUL", "SITA-SITA", "OSAT-HVUL", "TURA-TURA", "ZMAY-SBIC", "ZMAY-<7>", "<5>-HVUL", "BDIS-BDIS", "OSAT-OSAT", "TAES-OSAT", "<2>-OSAT", "LPER-LPER", "<3>-SITA", "<2>-<7>", "<1>-OSAT", "<4>-<8>", "SBIC-SITA", "SITA-HVUL", "LPER-HVUL", "<1>-<4>", "<1>-SITA", "BDIS-HVUL", "TAES-<8>", "TURA-OSAT", "OSAT-LPER", "OSAT-<5>", "<3>-LPER", "<8>-BDIS", "TURA-ATAU", "TURA-<1>", "<3>-<5>", "ATAU-ATAU", "<7>-BDIS", "SBIC-SBIC", "ATAU-TURA", "ATAU-<1>", "SBIC-BDIS", "TURA-SITA", "<3>-OSAT", "ATAU-LPER", "TURA-<5>", "<2>-ZMAY", "ZMAY-BDIS", "<4>-SITA", "<6>-SITA", "SBIC-ZMAY", "SBIC-<7>", "<8>-LPER", "<7>-SITA", "<7>-<8>", "OSAT-TURA", "<4>-LPER", "TAES-ZMAY", "<2>-SBIC", "OSAT-BDIS", "<8>-<5>", "<7>-TURA", "<5>-TURA", "<8>-TURA", "<7>-ATAU", "<5>-BDIS", "<8>-ATAU", "<3>-ZMAY", "<7>-LPER", "ZMAY-SITA", "HVUL-BDIS", "OSAT-ATAU", "SBIC-LPER", "TAES-SBIC", "SBIC-TURA", "LPER-BDIS", "<3>-SBIC", "<5>-ATAU", "SBIC-ATAU", "<7>-<7>", "TAES-<7>", "HVUL-TURA", "SITA-BDIS", "ZMAY-TURA", "ZMAY-LPER", "ZMAY-ATAU", "HVUL-ATAU", "LPER-OSAT", "LPER-<5>", "TAES-<6>", "SBIC-<8>", "SITA-LPER", "<4>-OSAT", "BDIS-LPER", "LPER-ATAU", "LPER-TURA", "SITA-TURA", "SITA-ATAU", "HVUL-LPER", "BDIS-ATAU", "ATAU-SITA", "BDIS-TURA", "ATAU-OSAT", "ZMAY-<8>", "<8>-OSAT", "<4>-ZMAY", "<3>-<3>", "<5>-SITA", "<3>-<7>", "<4>-SBIC", "<7>-OSAT", "<1>-<8>", "<4>-<6>", "<4>-<5>", "<6>-SBIC", "<6>-ZMAY", "OSAT-SITA", "<1>-ZMAY", "SBIC-OSAT", "TURA-ZMAY", "HVUL-SITA", "LPER-SITA", "<1>-SBIC", "ZMAY-OSAT", "<5>-<8>", "<5>-<5>", "TURA-SBIC", "TURA-<4>", "ATAU-<5>", "BDIS-SITA", "SITA-OSAT", "<4>-<7>", "<7>-<5>", "ATAU-<2>", "BDIS-OSAT", "HVUL-OSAT", "ATAU-<3>", "SBIC-<5>", "<6>-<7>", "ATAU-ZMAY", "<1>-<7>", "BDIS-<5>", "ZMAY-<5>", "<5>-ZMAY", "ATAU-SBIC", "TURA-<8>", "OSAT-ZMAY", "SITA-ZMAY", "SITA-SBIC", "<5>-SBIC", "<1>-<6>", "OSAT-SBIC", "HVUL-ZMAY", "LPER-ZMAY", "BDIS-ZMAY", "TURA-<7>", "HVUL-SBIC", "LPER-SBIC", "HVUL-<5>", "BDIS-SBIC", "SITA-<5>", "SITA-<7>", "SITA-<8>", "<8>-<4>", "<8>-<8>", "OSAT-<8>", "<5>-<7>", "ATAU-<7>", "ATAU-<4>", "ATAU-<8>", "OSAT-<7>", "LPER-<8>", "OSAT-<4>", "HVUL-<7>", "<7>-<3>", "LPER-<7>", "<8>-<3>", "HVUL-<8>", "BDIS-<7>", "<7>-<4>", "BDIS-<8>", "TURA-<6>", "<7>-<6>", "OSAT-<3>", "SBIC-<3>", "<4>-<4>", "ZMAY-<3>", "SBIC-<4>", "LPER-<4>", "HVUL-<2>", "HVUL-<3>", "BDIS-<4>", "BDIS-<3>", "<5>-<3>", "<5>-<6>", "<5>-<4>", "ZMAY-<4>", "OSAT-<6>", "LPER-<3>", "HVUL-<4>", "SITA-<3>", "SITA-<4>", "ZMAY-<6>", "SBIC-<6>", "<7>-<2>", "<8>-<2>", "LPER-<6>", "ATAU-<6>", "OSAT-<2>", "SBIC-<2>", "<5>-<2>", "ZMAY-<2>", "BDIS-<6>", "TURA-TAES", "<6>-<6>", "SITA-<6>", "TAES-TAES", "LPER-<2>", "BDIS-<2>", "SITA-<2>", "HVUL-<6>", "<1>-<1>", "<8>-<1>", "<7>-<1>", "OSAT-<1>", "HVUL-<1>", "<5>-<1>", "SBIC-<1>", "ZMAY-<1>", "LPER-<1>", "SITA-<1>", "BDIS-<1>", "<7>-TAES", "<8>-TAES", "HVUL-TAES", "OSAT-TAES", "SBIC-TAES", "ZMAY-TAES", "<5>-TAES", "LPER-TAES", "SITA-TAES", "BDIS-TAES"]}], {"autosize": false, "title": "GRAMPA Results: wheat_bd_B_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';

outfilename = "../../results/wheat/" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, p1=p1, p2=p2, footer=footer));