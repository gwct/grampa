############################################################
# For GRAMPA site, 12.19
# This generates the file "wheat_ab.html"
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

				<h1>Wheat Results  (AB sub-genomes)</h1>

				<div class="row">
					<a name="install"></a>
					<div class="col-1" id="jump_row">
						<div id="jump_container">
							<center>Files:
							<a class="jump_link" href="trees/wheat_A.tre" download>Species tree</a>
							<a class="jump_link" href="trees/wheat_AB_trees.txt" download>Gene trees</a>
							<a class="jump_link" href="output/wheat_AB_A_out.txt" download>Main output table</a>
							<a class="jump_link" href="wheat_all.html">All sub-genome results</a>
							<a class="jump_link" href="wheat_ad.html">AD sub-genome results</a>
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
						<h3>7,954 trees after filtering</h3>									
					</div>

					<div class="col-12-24">
						<center><h2>Optimal Wheat MUL-tree</h2></center>
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/wheat_AB_mul.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>	
						<h3>Optimal H1 node: TAES (A sub-genome)</h3>
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
pagefile = "wheat_ab.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Wheat (AB)"

head = RC.readResultsHead(title);
nav = RC.readNav(pagefile, "../../", "../yeast/", "");
footer = RC.readFooter();

p1 = '<div id="ffb57e6d-20e9-444e-ae61-e747556981bb" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("ffb57e6d-20e9-444e-ae61-e747556981bb", [{"opacity": 0.6, "y": [82588, 84298, 84298, 86060, 91215, 92545, 92545, 92559, 92820, 93285, 93558, 93633, 93865, 93908, 93924, 94215, 94215, 94601, 94691, 94811, 96225, 96225, 96278, 96354, 96585, 96624, 96665, 96704, 96907, 97112, 97327, 97329, 97404, 97484, 97484, 97501, 97634, 97939, 98021, 98197, 98211, 98235, 98253, 98352, 98380, 98392, 98411, 98450, 98482, 98570, 98579, 98582, 98646, 98688, 98690, 98698, 98775, 98783, 98823, 99136, 99314, 99343, 99409, 99433, 99436, 99461, 99489, 99491, 99491, 99510, 99536, 99537, 99659, 99669, 99686, 99686, 99699, 99713, 99717, 99718, 99742, 99772, 99777, 99788, 99818, 99828, 99849, 99859, 99861, 99897, 99899, 99907, 99928, 99929, 99966, 99984, 99996, 100001, 100001, 100028, 100031, 100052, 100059, 100086, 100104, 100132, 100168, 100198, 100266, 100270, 100279, 100285, 100290, 100299, 100299, 100303, 100303, 100309, 100317, 100343, 100364, 100364, 100377, 100381, 100388, 100394, 100406, 100414, 100423, 100434, 100447, 100450, 100454, 100460, 100469, 100483, 100483, 100487, 100490, 100495, 100507, 100520, 100521, 100521, 100535, 100546, 100583, 100591, 100597, 100624, 100625, 100625, 100635, 100635, 100651, 100691, 100697, 100701, 100732, 100767, 100768, 100778, 100784, 100807, 100839, 100847, 100854, 100859, 100862, 100897, 100925, 100948, 100948, 100951, 100953, 100963, 100965, 100995, 101006, 101022, 101037, 101054, 101069, 101077, 101089, 101131, 101144, 101149, 101177, 101272, 101317, 101317, 101323, 101329, 101393, 101433, 101444, 101494, 101507, 101536, 101540, 101552, 101572, 101576, 101580, 101590, 101610, 101643, 101659, 101675, 101683, 101693, 101705, 101777, 101780, 101807, 101823, 101843, 101875, 101885, 101885, 101982, 102239, 102248, 102281, 102407, 102425, 102450, 102474, 102507, 102512, 102592, 102615, 102616, 102622, 102643, 102658, 102665, 102770, 102801, 102860, 102918, 102929, 102956, 102992, 103014, 103015, 103015, 103047, 103052, 103065, 103065, 103097, 103097, 103126, 103158, 103301, 103333, 103459, 103459, 103469, 103604, 103618, 103636, 103764, 103789, 103838, 103840, 103869, 103937, 104047, 104054, 104098, 104160, 104164, 104204, 104280, 104380, 104432, 104438, 104461, 104512, 106889, 106889, 109075, 109192, 109265, 109269, 109318, 109366, 109417, 109472, 109627, 109668, 109682, 109697, 109724, 109741, 109817, 109892, 109965, 110205, 110284, 110301], "type": "scatter", "mode": "markers", "x": ["TAES-<2>", "<1>-ATAU", "<1>-<2>", "TAES-ATAU", "ST", "TAES-<1>", "TAES-TURA", "<9>-<9>", "<2>-<9>", "<3>-<9>", "<1>-<9>", "<1>-<3>", "<4>-<9>", "TAES-HVUL", "TAES-<3>", "<6>-<8>", "<6>-<9>", "<2>-<2>", "TAES-<9>", "<1>-HVUL", "<2>-<3>", "<2>-HVUL", "TURA-ATAU", "ZMAY-ZMAY", "<1>-<4>", "<2>-<4>", "ATAU-TURA", "TAES-BDIS", "<7>-<9>", "<1>-BDIS", "TURA-<9>", "ZMAY-<9>", "<5>-<9>", "<8>-<9>", "<8>-<6>", "TAES-<4>", "<1>-<5>", "<2>-<6>", "SBIC-<9>", "TAES-<5>", "<2>-<5>", "OSAT-<9>", "<2>-BDIS", "TAES-LPER", "<1>-<6>", "<1>-LPER", "<1>-SITA", "ATAU-<9>", "TURA-HVUL", "<2>-<8>", "TAES-SITA", "LPER-<9>", "HVUL-<9>", "<1>-<8>", "<1>-OSAT", "SITA-<9>", "BDIS-<9>", "<2>-SITA", "TAES-OSAT", "<3>-<6>", "TAES-<8>", "TURA-BDIS", "<8>-HVUL", "HVUL-HVUL", "<3>-<8>", "<7>-HVUL", "ATAU-HVUL", "<3>-BDIS", "<3>-<4>", "<2>-LPER", "TAES-<6>", "SBIC-HVUL", "SITA-SITA", "ZMAY-HVUL", "ZMAY-SBIC", "ZMAY-<7>", "BDIS-BDIS", "OSAT-HVUL", "<5>-HVUL", "<1>-<7>", "TURA-LPER", "TAES-ZMAY", "OSAT-OSAT", "<2>-OSAT", "LPER-LPER", "SBIC-SITA", "ATAU-ATAU", "<1>-ZMAY", "<3>-SITA", "<4>-<8>", "TAES-SBIC", "LPER-HVUL", "TURA-TURA", "SITA-HVUL", "TAES-<7>", "<1>-SBIC", "BDIS-HVUL", "OSAT-LPER", "OSAT-<5>", "<3>-LPER", "<2>-<7>", "<3>-<5>", "<8>-BDIS", "ATAU-BDIS", "SBIC-SBIC", "TURA-OSAT", "<7>-BDIS", "SBIC-BDIS", "<3>-OSAT", "<8>-ATAU", "TURA-SITA", "<7>-ATAU", "<6>-SITA", "<7>-SITA", "<7>-<8>", "SBIC-ZMAY", "SBIC-<7>", "<4>-SITA", "ZMAY-BDIS", "<8>-LPER", "OSAT-TURA", "SBIC-ATAU", "<4>-LPER", "<2>-ZMAY", "OSAT-ATAU", "<5>-ATAU", "OSAT-BDIS", "<5>-TURA", "<8>-<5>", "<5>-BDIS", "<7>-TURA", "ZMAY-SITA", "ZMAY-ATAU", "<8>-TURA", "TURA-<5>", "HVUL-BDIS", "<2>-SBIC", "HVUL-ATAU", "ATAU-LPER", "<7>-LPER", "HVUL-TURA", "<7>-<7>", "LPER-BDIS", "<3>-ZMAY", "SBIC-LPER", "SBIC-TURA", "<3>-SBIC", "TURA-<2>", "SITA-BDIS", "ZMAY-TURA", "LPER-ATAU", "ZMAY-LPER", "LPER-OSAT", "LPER-<5>", "SBIC-<8>", "TURA-<3>", "SITA-ATAU", "BDIS-ATAU", "SITA-LPER", "<4>-OSAT", "BDIS-LPER", "HVUL-LPER", "ZMAY-<8>", "LPER-TURA", "SITA-TURA", "<8>-OSAT", "<5>-SITA", "BDIS-TURA", "<3>-<3>", "<4>-ZMAY", "ATAU-SITA", "<4>-<6>", "<4>-<5>", "<4>-SBIC", "OSAT-SITA", "<7>-OSAT", "<3>-<7>", "<6>-SBIC", "<6>-ZMAY", "HVUL-SITA", "ATAU-OSAT", "SBIC-OSAT", "<5>-<8>", "TURA-ZMAY", "LPER-SITA", "ZMAY-OSAT", "BDIS-SITA", "<5>-<5>", "TURA-SBIC", "SITA-OSAT", "<4>-<7>", "HVUL-OSAT", "<7>-<5>", "BDIS-OSAT", "<6>-<7>", "SBIC-<5>", "ATAU-<5>", "<5>-ZMAY", "TURA-<4>", "OSAT-ZMAY", "BDIS-<5>", "ZMAY-<5>", "ATAU-ZMAY", "SITA-SBIC", "SITA-ZMAY", "<5>-SBIC", "HVUL-ZMAY", "OSAT-SBIC", "TURA-<8>", "ATAU-SBIC", "LPER-ZMAY", "HVUL-SBIC", "BDIS-ZMAY", "TURA-<7>", "LPER-SBIC", "BDIS-SBIC", "HVUL-<5>", "SITA-<5>", "<8>-<4>", "SITA-<7>", "SITA-<8>", "<8>-<8>", "ATAU-<3>", "OSAT-<8>", "<5>-<7>", "OSAT-<4>", "OSAT-<7>", "ATAU-<7>", "ATAU-<8>", "LPER-<8>", "HVUL-<7>", "<7>-<3>", "<8>-<3>", "LPER-<7>", "HVUL-<8>", "<7>-<4>", "BDIS-<7>", "BDIS-<8>", "<7>-<6>", "ATAU-<4>", "OSAT-<3>", "SBIC-<3>", "<4>-<4>", "TURA-<6>", "LPER-<4>", "SBIC-<4>", "HVUL-<2>", "HVUL-<3>", "ZMAY-<3>", "<5>-<3>", "BDIS-<4>", "BDIS-<3>", "<5>-<6>", "<5>-<4>", "OSAT-<6>", "ZMAY-<4>", "LPER-<3>", "HVUL-<4>", "ATAU-<1>", "ATAU-<2>", "SITA-<3>", "SITA-<4>", "ZMAY-<6>", "SBIC-<6>", "LPER-<6>", "<7>-<2>", "<8>-<2>", "ATAU-TAES", "<1>-<1>", "OSAT-<2>", "<5>-<2>", "SBIC-<2>", "ATAU-<6>", "BDIS-<6>", "ZMAY-<2>", "<6>-<6>", "SITA-<6>", "LPER-<2>", "BDIS-<2>", "HVUL-<6>", "TAES-TAES", "SITA-<2>", "TURA-TAES", "TURA-<1>", "HVUL-TAES", "<7>-TAES", "<8>-TAES", "OSAT-TAES", "<5>-TAES", "SBIC-TAES", "ZMAY-TAES", "HVUL-<1>", "LPER-TAES", "<7>-<1>", "BDIS-TAES", "SITA-TAES", "OSAT-<1>", "<8>-<1>", "<5>-<1>", "SBIC-<1>", "ZMAY-<1>", "LPER-<1>", "BDIS-<1>", "SITA-<1>"]}], {"autosize": false, "title": "GRAMPA Results: wheat_ab_A_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';
p2 = '<div id="3b7daf49-0dcd-4ef8-a084-1a3878b06025" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("3b7daf49-0dcd-4ef8-a084-1a3878b06025", [{"opacity": 0.6, "y": [84244, 94791, 95940, 95940, 98252, 98610, 99215, 100497, 100600, 100600, 100734, 100968, 101225, 101816, 102174, 102174, 103681, 104164, 104164, 104280, 104457, 104607, 104732, 104931, 105355, 105364, 105426, 105514, 105514, 105556, 105616, 105870, 106040, 106066, 106177, 106271, 106288, 106331, 106339, 106601, 106628, 106674, 106696, 106734, 106745, 106758, 106817, 106893, 107260, 107420, 107513, 107546, 107562, 107588, 107589, 107593, 107595, 107607, 107651, 107651, 107660, 107664, 107665, 107696, 107808, 107809, 107814, 107845, 107845, 107861, 107868, 107874, 107898, 107933, 107948, 107961, 107961, 107970, 107984, 107986, 108004, 108022, 108051, 108069, 108071, 108071, 108085, 108091, 108158, 108161, 108161, 108169, 108185, 108208, 108209, 108230, 108230, 108266, 108282, 108298, 108329, 108371, 108371, 108373, 108426, 108441, 108441, 108460, 108463, 108470, 108470, 108470, 108470, 108470, 108494, 108529, 108541, 108581, 108581, 108582, 108588, 108599, 108600, 108614, 108620, 108622, 108632, 108638, 108642, 108654, 108664, 108669, 108678, 108701, 108701, 108702, 108705, 108707, 108713, 108730, 108737, 108763, 108768, 108771, 108780, 108796, 108806, 108806, 108829, 108834, 108840, 108843, 108847, 108851, 108884, 108895, 108899, 108912, 108915, 108936, 108938, 108945, 108946, 108947, 109002, 109003, 109009, 109024, 109030, 109033, 109052, 109060, 109110, 109110, 109112, 109130, 109130, 109132, 109141, 109157, 109163, 109164, 109170, 109200, 109230, 109249, 109253, 109262, 109268, 109295, 109317, 109331, 109352, 109449, 109484, 109495, 109500, 109507, 109562, 109573, 109582, 109617, 109642, 109664, 109679, 109710, 109719, 109720, 109759, 109764, 109766, 109775, 109796, 109826, 109853, 109870, 109878, 109885, 109966, 109972, 109986, 109998, 110014, 110031, 110058, 110091, 110091, 110187, 110383, 110400, 110459, 110490, 110634, 110641, 110715, 110722, 110805, 110818, 110827, 110833, 110860, 110861, 110869, 110999, 111088, 111155, 111158, 111206, 111234, 111253, 111253, 111253, 111254, 111287, 111295, 111295, 111341, 111341, 111371, 111371, 111542, 111577, 111704, 111739, 111757, 111846, 111846, 111889, 112019, 112039, 112081, 112198, 112313, 112322, 112408, 112410, 112464, 112537, 112653, 112699, 112704, 112707, 112776, 115738, 116169, 116202, 116212, 116225, 116233, 116330, 116420, 116575, 116636, 116637, 117503, 117531, 117579, 117683, 117709, 117729, 117768, 118000, 118059, 118062], "type": "scatter", "mode": "markers", "x": ["TAES-TURA", "TAES-ATAU", "TAES-<1>", "TAES-<2>", "TURA-TAES", "ATAU-<2>", "ST", "<9>-<9>", "<1>-<2>", "<1>-TAES", "<2>-<9>", "<2>-<2>", "<3>-<9>", "<4>-<9>", "<6>-<8>", "<6>-<9>", "TAES-<9>", "<2>-<3>", "<2>-HVUL", "<1>-<9>", "ZMAY-ZMAY", "TAES-HVUL", "<2>-<4>", "<7>-<9>", "ZMAY-<9>", "TURA-<9>", "<5>-<9>", "<8>-<9>", "<8>-<6>", "ATAU-HVUL", "<1>-HVUL", "TAES-<3>", "<2>-<6>", "SBIC-<9>", "ATAU-<9>", "OSAT-<9>", "<1>-<3>", "<2>-<5>", "<2>-BDIS", "TURA-HVUL", "LPER-<9>", "<2>-<8>", "HVUL-<9>", "<1>-BDIS", "SITA-<9>", "TAES-BDIS", "BDIS-<9>", "<2>-SITA", "<3>-<6>", "ATAU-BDIS", "TURA-BDIS", "<8>-HVUL", "<3>-<8>", "TAES-SITA", "TAES-LPER", "HVUL-HVUL", "<1>-LPER", "<7>-HVUL", "<3>-BDIS", "<3>-<4>", "<1>-<5>", "TAES-<5>", "<2>-LPER", "SBIC-HVUL", "ZMAY-HVUL", "SITA-SITA", "TAES-<4>", "ZMAY-SBIC", "ZMAY-<7>", "BDIS-BDIS", "OSAT-HVUL", "<5>-HVUL", "TURA-LPER", "OSAT-OSAT", "<2>-OSAT", "ATAU-<3>", "<1>-<4>", "LPER-LPER", "SBIC-SITA", "<1>-OSAT", "<3>-SITA", "<4>-<8>", "ATAU-ATAU", "LPER-HVUL", "<1>-SITA", "TAES-OSAT", "SITA-HVUL", "TURA-TURA", "BDIS-HVUL", "OSAT-LPER", "OSAT-<5>", "<2>-<7>", "<3>-LPER", "<3>-<5>", "<8>-BDIS", "TURA-ATAU", "TURA-<1>", "SBIC-SBIC", "ATAU-LPER", "TURA-OSAT", "<7>-BDIS", "ATAU-TURA", "ATAU-<1>", "SBIC-BDIS", "<3>-OSAT", "TURA-SITA", "<6>-SITA", "<4>-SITA", "TAES-<8>", "ZMAY-BDIS", "SBIC-ZMAY", "<7>-SITA", "<7>-<8>", "SBIC-<7>", "<8>-LPER", "<2>-ZMAY", "<4>-LPER", "OSAT-TURA", "OSAT-BDIS", "<8>-<5>", "OSAT-ATAU", "<5>-BDIS", "<7>-ATAU", "<8>-ATAU", "ZMAY-SITA", "<5>-ATAU", "<2>-SBIC", "TURA-<5>", "<5>-TURA", "<7>-LPER", "HVUL-BDIS", "<3>-ZMAY", "SBIC-ATAU", "<7>-TURA", "LPER-BDIS", "SBIC-LPER", "TAES-ZMAY", "<7>-<7>", "HVUL-ATAU", "<8>-TURA", "<3>-SBIC", "ZMAY-ATAU", "TURA-<2>", "SITA-BDIS", "ZMAY-LPER", "SBIC-TURA", "LPER-OSAT", "LPER-<5>", "SBIC-<8>", "TAES-SBIC", "ATAU-OSAT", "ATAU-SITA", "HVUL-TURA", "LPER-ATAU", "ZMAY-TURA", "SITA-ATAU", "SITA-LPER", "BDIS-ATAU", "<3>-<3>", "BDIS-LPER", "<4>-OSAT", "TURA-<3>", "HVUL-LPER", "ZMAY-<8>", "LPER-TURA", "TAES-<7>", "<8>-OSAT", "SITA-TURA", "<1>-<8>", "<5>-SITA", "<4>-ZMAY", "BDIS-TURA", "<4>-<6>", "<4>-<5>", "<4>-SBIC", "<3>-<7>", "<7>-OSAT", "OSAT-SITA", "ATAU-<5>", "<1>-ZMAY", "<6>-SBIC", "TAES-<6>", "<6>-ZMAY", "HVUL-SITA", "SBIC-OSAT", "<1>-SBIC", "TURA-ZMAY", "<5>-<8>", "LPER-SITA", "ZMAY-OSAT", "BDIS-SITA", "<5>-<5>", "TURA-SBIC", "SITA-OSAT", "<4>-<7>", "HVUL-OSAT", "<7>-<5>", "BDIS-OSAT", "ATAU-ZMAY", "<6>-<7>", "<1>-<7>", "SBIC-<5>", "<1>-<6>", "ATAU-SBIC", "<5>-ZMAY", "BDIS-<5>", "ZMAY-<5>", "OSAT-ZMAY", "SITA-SBIC", "TURA-<4>", "SITA-ZMAY", "<5>-SBIC", "HVUL-ZMAY", "OSAT-SBIC", "TURA-<8>", "LPER-ZMAY", "HVUL-SBIC", "BDIS-ZMAY", "LPER-SBIC", "TURA-<7>", "BDIS-SBIC", "ATAU-<4>", "HVUL-<5>", "SITA-<5>", "<8>-<4>", "SITA-<7>", "SITA-<8>", "<8>-<8>", "ATAU-<8>", "ATAU-<7>", "OSAT-<8>", "<5>-<7>", "OSAT-<7>", "OSAT-<4>", "LPER-<8>", "HVUL-<7>", "<7>-<3>", "<8>-<3>", "LPER-<7>", "HVUL-<8>", "<7>-<4>", "BDIS-<7>", "BDIS-<8>", "<7>-<6>", "OSAT-<3>", "SBIC-<3>", "<4>-<4>", "TURA-<6>", "LPER-<4>", "HVUL-<2>", "ZMAY-<3>", "HVUL-<3>", "SBIC-<4>", "<5>-<3>", "BDIS-<4>", "BDIS-<3>", "<5>-<6>", "<5>-<4>", "ZMAY-<4>", "OSAT-<6>", "LPER-<3>", "HVUL-<4>", "SITA-<3>", "ATAU-<6>", "ATAU-TAES", "SITA-<4>", "ZMAY-<6>", "SBIC-<6>", "LPER-<6>", "<7>-<2>", "<8>-<2>", "OSAT-<2>", "<5>-<2>", "SBIC-<2>", "BDIS-<6>", "ZMAY-<2>", "<6>-<6>", "SITA-<6>", "LPER-<2>", "BDIS-<2>", "TAES-TAES", "HVUL-<6>", "SITA-<2>", "<1>-<1>", "OSAT-<1>", "<8>-<1>", "<7>-<1>", "HVUL-<1>", "<5>-<1>", "SBIC-<1>", "ZMAY-<1>", "LPER-<1>", "SITA-<1>", "BDIS-<1>", "HVUL-TAES", "<7>-TAES", "<8>-TAES", "OSAT-TAES", "SBIC-TAES", "<5>-TAES", "ZMAY-TAES", "LPER-TAES", "BDIS-TAES", "SITA-TAES"]}], {"autosize": false, "title": "GRAMPA Results: wheat_ab_B_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';

outfilename = "../../results/wheat/" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, p1=p1, p2=p2, footer=footer));