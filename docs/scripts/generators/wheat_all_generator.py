############################################################
# For GRAMPA site, 12.19
# This generates the file "wheat_all.html"
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

				<h1>Wheat Results (all sub-genomes)</h1>

				<div class="row">
					<a name="install"></a>
					<div class="col-1" id="jump_row">
						<div id="jump_container">
							<center>Files:
							<a class="jump_link" href="trees/wheat_A.tre" download>Species tree</a>
							<a class="jump_link" href="trees/wheat_all_trees.txt" download>Gene trees</a>
							<a class="jump_link" href="output/wheat_all_A_out.txt" download>Main output table</a>
							<a class="jump_link" href="wheat_ab.html">AB sub-genome results</a>
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
						<h3>7,362 trees after filtering</h3>									
					</div>

					<div class="col-12-24">
						<center><h2>Optimal Wheat MUL-tree</h2></center>
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/wheat_all_mul.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>	
						<h3>Optimal H1 node: TAES (A sub-genome)</h3>
						<h3>Optimal H2 node: 2 (B sub-genome)</h3>
						<h3>Second H2 node: ATAU (D sub-genome)</h3>				
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

				<div class="row">
					<div class="col-1">
						<h2>Score plot with alternate species topology 2</h2>
						<center>
						{p3}
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
pagefile = "wheat_all.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Wheat (all)"

head = RC.readResultsHead(title);
nav = RC.readNav(pagefile, "../../", "../yeast/", "");
footer = RC.readFooter();

p1 = '<div id="ce2743d9-bab8-432f-a19f-7568c2b7f136" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("ce2743d9-bab8-432f-a19f-7568c2b7f136", [{"opacity": 0.6, "y": [75321, 78413, 82751, 82751, 92517, 92517, 92947, 93161, 93394, 93668, 93868, 94395, 94753, 94948, 95042, 95180, 95216, 95531, 95606, 96077, 96077, 96628, 96783, 96957, 96957, 97055, 97336, 97340, 97421, 97746, 98034, 98035, 98157, 98222, 98451, 98560, 98573, 98614, 98650, 98663, 98687, 98792, 98816, 98868, 98888, 98892, 98934, 98934, 98961, 99012, 99033, 99093, 99093, 99186, 99445, 99473, 99707, 99708, 99746, 99759, 99763, 99806, 99809, 99870, 99870, 99892, 99902, 99952, 99952, 99966, 99968, 100013, 100015, 100021, 100039, 100051, 100053, 100061, 100064, 100067, 100108, 100132, 100142, 100148, 100151, 100154, 100154, 100156, 100160, 100207, 100216, 100225, 100228, 100240, 100244, 100263, 100286, 100319, 100326, 100377, 100386, 100394, 100394, 100430, 100435, 100446, 100463, 100476, 100477, 100487, 100503, 100504, 100532, 100534, 100542, 100549, 100562, 100571, 100585, 100604, 100604, 100607, 100607, 100611, 100611, 100620, 100624, 100626, 100644, 100644, 100679, 100699, 100700, 100706, 100714, 100727, 100729, 100731, 100747, 100754, 100769, 100775, 100779, 100785, 100790, 100811, 100828, 100832, 100834, 100853, 100862, 100868, 100872, 100895, 100915, 100928, 100928, 100934, 100939, 100943, 100962, 100994, 101020, 101022, 101039, 101054, 101064, 101097, 101098, 101117, 101140, 101153, 101165, 101165, 101177, 101180, 101207, 101215, 101239, 101244, 101272, 101280, 101282, 101297, 101330, 101332, 101358, 101364, 101388, 101396, 101451, 101467, 101506, 101510, 101532, 101547, 101626, 101630, 101673, 101676, 101715, 101719, 101738, 101741, 101747, 101754, 101754, 101775, 101787, 101792, 101793, 101794, 101859, 101864, 101876, 101889, 101910, 101947, 101947, 101966, 101984, 102055, 102242, 102304, 102308, 102428, 102498, 102529, 102538, 102545, 102569, 102607, 102615, 102639, 102663, 102670, 102726, 102767, 102790, 102798, 102906, 102979, 103025, 103036, 103057, 103084, 103096, 103099, 103126, 103131, 103136, 103136, 103140, 103140, 103180, 103180, 103227, 103324, 103430, 103468, 103552, 103555, 103574, 103617, 103662, 103662, 103719, 103804, 103805, 103916, 103962, 103988, 103998, 104053, 104056, 104126, 104238, 104275, 104296, 104355, 107054, 107255, 110905, 110905, 112961, 113043, 113063, 113144, 113152, 113287, 113325, 113386, 113418, 113476, 113497, 113525, 113537, 113613, 113625, 113630, 113647, 113842, 113905, 113906], "type": "scatter", "mode": "markers", "x": ["TAES-ATAU", "TAES-<2>", "<1>-ATAU", "<1>-<2>", "TAES-<1>", "TAES-TURA", "TAES-<3>", "ST", "<2>-<2>", "<1>-<3>", "TAES-HVUL", "<2>-<9>", "<1>-<9>", "<3>-<9>", "<9>-<9>", "ATAU-TURA", "<1>-HVUL", "TAES-<9>", "<4>-<9>", "<6>-<8>", "<6>-<9>", "<1>-<4>", "TAES-BDIS", "<2>-<3>", "<2>-HVUL", "TAES-<4>", "<2>-<4>", "ZMAY-ZMAY", "<1>-BDIS", "<1>-<5>", "TURA-ATAU", "TAES-<5>", "<1>-<6>", "<2>-<6>", "TAES-LPER", "<7>-<9>", "<2>-<5>", "TAES-SITA", "<1>-SITA", "<1>-LPER", "<1>-<8>", "<2>-BDIS", "TAES-OSAT", "<1>-OSAT", "<2>-<8>", "TAES-<6>", "TURA-HVUL", "TURA-<9>", "ZMAY-<9>", "<5>-<9>", "TAES-<8>", "<8>-<9>", "<8>-<6>", "<2>-SITA", "<3>-<6>", "SBIC-<9>", "OSAT-<9>", "<1>-<7>", "TAES-<7>", "TURA-BDIS", "TAES-ZMAY", "<3>-<8>", "HVUL-HVUL", "<8>-HVUL", "<2>-LPER", "<7>-HVUL", "TAES-SBIC", "<3>-BDIS", "<3>-<4>", "LPER-<9>", "SBIC-HVUL", "<1>-ZMAY", "ATAU-<9>", "ATAU-HVUL", "SITA-SITA", "ZMAY-HVUL", "OSAT-HVUL", "<2>-OSAT", "<5>-HVUL", "SITA-<9>", "HVUL-<9>", "BDIS-BDIS", "TURA-LPER", "<1>-SBIC", "BDIS-<9>", "ZMAY-SBIC", "ZMAY-<7>", "OSAT-OSAT", "TURA-TURA", "<2>-<7>", "LPER-LPER", "ATAU-ATAU", "SBIC-SITA", "<3>-SITA", "LPER-HVUL", "SITA-HVUL", "ATAU-TAES", "BDIS-HVUL", "<3>-<5>", "<3>-LPER", "<4>-<8>", "OSAT-LPER", "OSAT-<5>", "<8>-BDIS", "TURA-OSAT", "TURA-<2>", "SBIC-SBIC", "OSAT-TURA", "<7>-ATAU", "<8>-ATAU", "<7>-BDIS", "ATAU-BDIS", "<5>-TURA", "<3>-OSAT", "SBIC-BDIS", "OSAT-ATAU", "<5>-ATAU", "SBIC-ATAU", "<7>-TURA", "SBIC-ZMAY", "SBIC-<7>", "<6>-SITA", "<2>-ZMAY", "<8>-TURA", "TURA-SITA", "<4>-SITA", "ZMAY-ATAU", "ZMAY-BDIS", "<7>-SITA", "<7>-<8>", "OSAT-BDIS", "<4>-LPER", "HVUL-TURA", "SBIC-TURA", "<5>-BDIS", "<8>-LPER", "ZMAY-SITA", "<2>-SBIC", "ZMAY-TURA", "HVUL-ATAU", "LPER-ATAU", "<3>-ZMAY", "TURA-<5>", "LPER-BDIS", "HVUL-BDIS", "<7>-<7>", "<7>-LPER", "BDIS-ATAU", "SITA-ATAU", "SITA-BDIS", "<3>-SBIC", "SBIC-LPER", "ATAU-LPER", "<8>-<5>", "LPER-TURA", "LPER-OSAT", "LPER-<5>", "SBIC-<8>", "ZMAY-LPER", "SITA-TURA", "BDIS-TURA", "<4>-OSAT", "<3>-<3>", "SITA-LPER", "BDIS-LPER", "TURA-<3>", "HVUL-LPER", "<5>-SITA", "<4>-ZMAY", "ZMAY-<8>", "<8>-OSAT", "<3>-<7>", "<4>-<6>", "<4>-<5>", "<4>-SBIC", "OSAT-SITA", "<6>-ZMAY", "<7>-OSAT", "ATAU-SITA", "<6>-SBIC", "HVUL-SITA", "LPER-SITA", "TURA-ZMAY", "SBIC-OSAT", "ATAU-OSAT", "BDIS-SITA", "ZMAY-OSAT", "<5>-<5>", "<5>-<8>", "TURA-SBIC", "<4>-<7>", "SITA-OSAT", "BDIS-OSAT", "HVUL-OSAT", "<6>-<7>", "<7>-<5>", "<5>-ZMAY", "SBIC-<5>", "OSAT-ZMAY", "SITA-ZMAY", "BDIS-<5>", "SITA-SBIC", "<5>-SBIC", "ZMAY-<5>", "ATAU-ZMAY", "HVUL-ZMAY", "ATAU-<5>", "LPER-ZMAY", "TURA-<8>", "TURA-<4>", "BDIS-ZMAY", "OSAT-SBIC", "HVUL-SBIC", "ATAU-SBIC", "TURA-<7>", "LPER-SBIC", "BDIS-SBIC", "SITA-<7>", "SITA-<8>", "SITA-<5>", "HVUL-<5>", "<8>-<8>", "<8>-<4>", "<5>-<7>", "OSAT-<8>", "OSAT-<7>", "LPER-<8>", "ATAU-<7>", "HVUL-<7>", "OSAT-<4>", "LPER-<7>", "BDIS-<7>", "ATAU-<8>", "BDIS-<8>", "HVUL-<8>", "ATAU-<3>", "<7>-<3>", "<7>-<4>", "<8>-<3>", "<7>-<6>", "OSAT-<3>", "TURA-<6>", "SBIC-<3>", "<4>-<4>", "LPER-<4>", "<5>-<3>", "SBIC-<4>", "ATAU-<4>", "ZMAY-<3>", "OSAT-<6>", "<5>-<6>", "<5>-<4>", "BDIS-<4>", "BDIS-<3>", "HVUL-<2>", "HVUL-<3>", "ZMAY-<4>", "LPER-<3>", "HVUL-<4>", "SITA-<3>", "ZMAY-<6>", "SBIC-<6>", "SITA-<4>", "LPER-<6>", "ATAU-<1>", "ATAU-<2>", "<7>-<2>", "OSAT-<2>", "<8>-<2>", "<5>-<2>", "<6>-<6>", "SBIC-<2>", "BDIS-<6>", "SITA-<6>", "ZMAY-<2>", "ATAU-<6>", "LPER-<2>", "BDIS-<2>", "HVUL-<6>", "SITA-<2>", "<1>-<1>", "TAES-TAES", "TURA-TAES", "TURA-<1>", "HVUL-<1>", "<7>-<1>", "OSAT-<1>", "<8>-<1>", "<5>-<1>", "SBIC-<1>", "ZMAY-<1>", "HVUL-TAES", "<7>-TAES", "OSAT-TAES", "<8>-TAES", "<5>-TAES", "LPER-<1>", "BDIS-<1>", "SITA-<1>", "SBIC-TAES", "ZMAY-TAES", "LPER-TAES", "SITA-TAES", "BDIS-TAES"]}], {"autosize": false, "title": "GRAMPA Results: wheat_all_A_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>'
p2 = '<div id="69b36e1d-7c73-49d3-b0cb-32aafccc34e0" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("69b36e1d-7c73-49d3-b0cb-32aafccc34e0", [{"opacity": 0.6, "y": [79990, 80558, 82435, 82435, 87899, 89739, 90921, 91483, 91543, 92151, 92636, 92636, 93189, 93249, 93462, 93462, 94022, 94030, 94046, 94886, 94978, 95187, 95330, 95468, 95514, 95572, 95588, 95593, 95622, 95627, 95656, 95728, 95728, 95932, 95935, 96080, 96111, 96112, 96168, 96192, 96192, 96207, 96294, 96338, 96364, 96485, 96537, 96558, 96578, 96606, 96614, 96627, 96644, 96654, 96705, 96705, 96716, 96726, 96741, 96743, 96769, 96784, 96787, 96792, 96820, 96823, 96833, 96845, 96852, 96888, 96914, 96917, 96917, 96920, 96955, 96955, 96980, 96984, 96986, 97008, 97012, 97027, 97036, 97066, 97074, 97092, 97104, 97104, 97108, 97140, 97148, 97164, 97164, 97197, 97199, 97220, 97234, 97235, 97263, 97263, 97275, 97308, 97315, 97340, 97364, 97364, 97365, 97374, 97383, 97383, 97390, 97394, 97403, 97403, 97403, 97408, 97417, 97422, 97437, 97437, 97455, 97459, 97472, 97489, 97492, 97493, 97494, 97506, 97521, 97522, 97525, 97527, 97549, 97561, 97561, 97568, 97574, 97578, 97587, 97590, 97601, 97608, 97620, 97626, 97636, 97652, 97653, 97660, 97662, 97669, 97670, 97687, 97705, 97705, 97720, 97721, 97740, 97762, 97777, 97786, 97802, 97820, 97847, 97878, 97878, 97880, 97883, 97900, 97903, 97915, 97946, 97951, 97951, 97965, 97972, 97984, 97993, 98030, 98057, 98063, 98067, 98071, 98078, 98114, 98129, 98145, 98146, 98171, 98199, 98199, 98250, 98254, 98272, 98295, 98300, 98315, 98332, 98396, 98407, 98410, 98425, 98458, 98463, 98506, 98518, 98524, 98535, 98536, 98545, 98560, 98562, 98578, 98593, 98593, 98632, 98668, 98685, 98693, 98713, 98744, 98744, 98762, 98783, 98844, 98868, 99037, 99112, 99114, 99198, 99239, 99263, 99297, 99318, 99353, 99364, 99381, 99418, 99444, 99474, 99535, 99589, 99595, 99620, 99726, 99830, 99842, 99856, 99881, 99906, 99911, 99947, 99957, 99957, 99965, 99965, 99965, 100003, 100003, 100037, 100059, 100144, 100262, 100285, 100386, 100402, 100405, 100459, 100567, 100636, 100649, 100662, 100709, 100775, 100823, 100834, 100839, 100909, 100913, 101094, 101136, 101156, 101208, 103577, 103859, 103918, 103932, 103942, 104036, 104062, 104121, 104182, 104228, 104282, 104291, 110532, 110587, 110590, 110648, 110695, 110719, 110773, 110957, 111015, 111021], "type": "scatter", "mode": "markers", "x": ["TAES-ATAU", "TAES-TURA", "TAES-<1>", "TAES-<2>", "<2>-<2>", "ST", "<2>-<9>", "<3>-<9>", "<9>-<9>", "<4>-<9>", "<6>-<8>", "<6>-<9>", "TAES-<9>", "TAES-HVUL", "<2>-<3>", "<2>-HVUL", "ZMAY-ZMAY", "TAES-<3>", "<2>-<4>", "<1>-<9>", "<2>-<6>", "<7>-<9>", "<2>-<5>", "<2>-BDIS", "TAES-BDIS", "<1>-HVUL", "TURA-<9>", "ZMAY-<9>", "<5>-<9>", "<2>-<8>", "TURA-HVUL", "<8>-<9>", "<8>-<6>", "<2>-SITA", "ATAU-HVUL", "TAES-<4>", "<1>-BDIS", "SBIC-<9>", "TAES-<5>", "<1>-<2>", "<1>-TAES", "<3>-<6>", "TAES-SITA", "OSAT-<9>", "TAES-LPER", "ATAU-<9>", "TURA-BDIS", "<3>-<8>", "HVUL-HVUL", "LPER-<9>", "<1>-<3>", "<8>-HVUL", "<2>-LPER", "<7>-HVUL", "<3>-BDIS", "<3>-<4>", "SITA-<9>", "SBIC-HVUL", "TAES-OSAT", "HVUL-<9>", "<1>-LPER", "BDIS-<9>", "ATAU-BDIS", "SITA-SITA", "ZMAY-HVUL", "OSAT-HVUL", "<5>-HVUL", "<2>-OSAT", "TAES-<8>", "BDIS-BDIS", "OSAT-OSAT", "ZMAY-SBIC", "ZMAY-<7>", "TURA-LPER", "<1>-<5>", "TURA-TURA", "<2>-<7>", "SBIC-SITA", "LPER-LPER", "<3>-SITA", "LPER-HVUL", "SITA-HVUL", "ATAU-ATAU", "TURA-<2>", "<1>-OSAT", "BDIS-HVUL", "TURA-ATAU", "TURA-<1>", "<3>-<5>", "<4>-<8>", "<3>-LPER", "OSAT-LPER", "OSAT-<5>", "<8>-BDIS", "TAES-<6>", "TURA-OSAT", "SBIC-SBIC", "<1>-SITA", "ATAU-TURA", "ATAU-<1>", "<7>-BDIS", "SBIC-BDIS", "<3>-OSAT", "OSAT-TURA", "SBIC-ZMAY", "SBIC-<7>", "TAES-ZMAY", "<6>-SITA", "OSAT-ATAU", "<2>-ZMAY", "TURA-SITA", "<4>-SITA", "ZMAY-BDIS", "<7>-SITA", "<7>-<8>", "<5>-TURA", "<7>-ATAU", "<5>-ATAU", "<8>-ATAU", "ATAU-LPER", "TAES-<7>", "OSAT-BDIS", "<4>-LPER", "SBIC-ATAU", "<7>-TURA", "<8>-LPER", "<5>-BDIS", "ZMAY-SITA", "TAES-SBIC", "<8>-TURA", "<2>-SBIC", "<1>-<4>", "<3>-ZMAY", "ZMAY-ATAU", "LPER-BDIS", "TURA-<5>", "HVUL-BDIS", "SBIC-TURA", "<7>-<7>", "HVUL-ATAU", "<7>-LPER", "<3>-<3>", "LPER-ATAU", "SITA-BDIS", "SBIC-LPER", "<3>-SBIC", "ZMAY-TURA", "<8>-<5>", "SITA-ATAU", "HVUL-TURA", "BDIS-ATAU", "SBIC-<8>", "LPER-OSAT", "LPER-<5>", "ATAU-<2>", "ZMAY-LPER", "LPER-TURA", "SITA-TURA", "<4>-OSAT", "BDIS-TURA", "SITA-LPER", "BDIS-LPER", "HVUL-LPER", "ATAU-SITA", "<4>-ZMAY", "<5>-SITA", "TURA-<3>", "ZMAY-<8>", "ATAU-OSAT", "<8>-OSAT", "<3>-<7>", "<4>-<6>", "<4>-<5>", "OSAT-SITA", "<4>-SBIC", "<6>-ZMAY", "<7>-OSAT", "<6>-SBIC", "HVUL-SITA", "LPER-SITA", "TURA-ZMAY", "SBIC-OSAT", "<1>-ZMAY", "BDIS-SITA", "<1>-<8>", "<5>-<5>", "ZMAY-OSAT", "<5>-<8>", "TURA-SBIC", "<1>-SBIC", "<4>-<7>", "SITA-OSAT", "ATAU-<5>", "BDIS-OSAT", "HVUL-OSAT", "<6>-<7>", "<7>-<5>", "ATAU-<3>", "SBIC-<5>", "<5>-ZMAY", "ATAU-ZMAY", "OSAT-ZMAY", "SITA-ZMAY", "BDIS-<5>", "<1>-<7>", "SITA-SBIC", "<5>-SBIC", "ZMAY-<5>", "HVUL-ZMAY", "ATAU-SBIC", "LPER-ZMAY", "BDIS-ZMAY", "OSAT-SBIC", "TURA-<8>", "TURA-<4>", "HVUL-SBIC", "TURA-<7>", "LPER-SBIC", "BDIS-SBIC", "SITA-<7>", "SITA-<8>", "SITA-<5>", "HVUL-<5>", "<8>-<8>", "<1>-<6>", "<8>-<4>", "<5>-<7>", "OSAT-<8>", "ATAU-<7>", "OSAT-<7>", "ATAU-<8>", "LPER-<8>", "ATAU-<4>", "HVUL-<7>", "OSAT-<4>", "LPER-<7>", "BDIS-<7>", "BDIS-<8>", "HVUL-<8>", "<7>-<3>", "<7>-<4>", "<8>-<3>", "<7>-<6>", "OSAT-<3>", "SBIC-<3>", "TURA-<6>", "<4>-<4>", "LPER-<4>", "<5>-<3>", "SBIC-<4>", "ZMAY-<3>", "BDIS-<4>", "BDIS-<3>", "<5>-<6>", "OSAT-<6>", "<5>-<4>", "HVUL-<2>", "HVUL-<3>", "ATAU-TAES", "ZMAY-<4>", "LPER-<3>", "HVUL-<4>", "SITA-<3>", "SBIC-<6>", "SITA-<4>", "ZMAY-<6>", "LPER-<6>", "<7>-<2>", "ATAU-<6>", "<8>-<2>", "OSAT-<2>", "TURA-TAES", "<5>-<2>", "<6>-<6>", "SBIC-<2>", "BDIS-<6>", "SITA-<6>", "ZMAY-<2>", "LPER-<2>", "BDIS-<2>", "HVUL-<6>", "SITA-<2>", "<1>-<1>", "OSAT-<1>", "<5>-<1>", "<7>-<1>", "<8>-<1>", "SBIC-<1>", "HVUL-<1>", "ZMAY-<1>", "TAES-TAES", "LPER-<1>", "BDIS-<1>", "SITA-<1>", "<7>-TAES", "HVUL-TAES", "<8>-TAES", "OSAT-TAES", "<5>-TAES", "SBIC-TAES", "ZMAY-TAES", "LPER-TAES", "SITA-TAES", "BDIS-TAES"]}], {"autosize": false, "title": "GRAMPA Results: wheat_all_B_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';
p3 = '<div id="b92a21d3-745b-422e-8f0d-b76efdacdfa4" style="height: 500; width: 1000px;" class="plotly-graph-div"></div><script type="text/javascript">window.PLOTLYENV=window.PLOTLYENV || {};window.PLOTLYENV.BASE_URL="https://plot.ly";Plotly.newPlot("b92a21d3-745b-422e-8f0d-b76efdacdfa4", [{"opacity": 0.6, "y": [75432, 77978, 82455, 82455, 90861, 90861, 92307, 92446, 92665, 93236, 93323, 93659, 94203, 94339, 94551, 94587, 94720, 94821, 94875, 95366, 95366, 96094, 96244, 96244, 96315, 96563, 96600, 96636, 96877, 97325, 97448, 97522, 97618, 97760, 97849, 97871, 97896, 98104, 98108, 98125, 98169, 98174, 98254, 98256, 98299, 98319, 98352, 98381, 98381, 98437, 98471, 98488, 98605, 98722, 98732, 98756, 98999, 99024, 99033, 99086, 99109, 99128, 99164, 99167, 99183, 99186, 99231, 99231, 99245, 99262, 99280, 99311, 99326, 99353, 99353, 99355, 99357, 99360, 99366, 99383, 99404, 99426, 99439, 99439, 99439, 99441, 99490, 99504, 99511, 99522, 99525, 99527, 99544, 99562, 99616, 99633, 99663, 99673, 99682, 99682, 99682, 99708, 99745, 99746, 99780, 99812, 99814, 99815, 99826, 99835, 99849, 99854, 99856, 99867, 99879, 99879, 99885, 99887, 99903, 99908, 99914, 99917, 99925, 99933, 99933, 99954, 99955, 99958, 99961, 99971, 99982, 100009, 100015, 100016, 100024, 100026, 100046, 100063, 100079, 100081, 100095, 100106, 100117, 100130, 100136, 100143, 100149, 100151, 100156, 100182, 100196, 100209, 100209, 100213, 100217, 100225, 100229, 100248, 100279, 100297, 100300, 100310, 100331, 100358, 100369, 100382, 100387, 100405, 100415, 100428, 100436, 100452, 100452, 100459, 100466, 100475, 100501, 100524, 100548, 100554, 100555, 100565, 100583, 100611, 100633, 100645, 100681, 100736, 100754, 100777, 100795, 100802, 100805, 100814, 100827, 100832, 100899, 100906, 100917, 100939, 100945, 100952, 101010, 101010, 101019, 101024, 101028, 101046, 101049, 101063, 101079, 101145, 101175, 101191, 101230, 101230, 101242, 101265, 101342, 101440, 101480, 101543, 101586, 101595, 101674, 101709, 101749, 101771, 101774, 101782, 101817, 101837, 101851, 101883, 101915, 101935, 101941, 102002, 102059, 102068, 102082, 102192, 102296, 102320, 102346, 102367, 102379, 102410, 102422, 102422, 102428, 102428, 102428, 102468, 102468, 102515, 102602, 102716, 102746, 102811, 102811, 102836, 102845, 102854, 102855, 102904, 103008, 103072, 103098, 103104, 103213, 103252, 103276, 103277, 103336, 103352, 103530, 103580, 103583, 103648, 105276, 106523, 109918, 109918, 111700, 111851, 111917, 111985, 112029, 112073, 112124, 112290, 112355, 112384, 112696, 112705, 112778, 112861, 112883, 112909, 112938, 113109, 113171, 113183], "type": "scatter", "mode": "markers", "x": ["TAES-TURA", "TAES-<2>", "<1>-<2>", "<1>-TURA", "TAES-ATAU", "TAES-<1>", "TAES-<3>", "ST", "<2>-<2>", "TAES-HVUL", "<1>-<3>", "<2>-<9>", "<3>-<9>", "<9>-<9>", "<1>-HVUL", "TURA-ATAU", "<1>-<9>", "TAES-<9>", "<4>-<9>", "<6>-<8>", "<6>-<9>", "TAES-BDIS", "<2>-<3>", "<2>-HVUL", "TAES-<4>", "<1>-<4>", "ZMAY-ZMAY", "<2>-<4>", "<1>-BDIS", "TAES-<5>", "ATAU-TURA", "<2>-<6>", "<1>-<5>", "TAES-LPER", "<7>-<9>", "<2>-<5>", "TAES-SITA", "<2>-BDIS", "<1>-<6>", "TAES-OSAT", "TAES-<6>", "<2>-<8>", "ZMAY-<9>", "<1>-SITA", "<5>-<9>", "TAES-<8>", "<1>-LPER", "<8>-<9>", "<8>-<6>", "<1>-<8>", "<2>-SITA", "ATAU-HVUL", "TURA-<9>", "<1>-OSAT", "<3>-<6>", "SBIC-<9>", "OSAT-<9>", "TAES-<7>", "TAES-ZMAY", "<3>-<8>", "HVUL-HVUL", "ATAU-<9>", "<2>-LPER", "<8>-HVUL", "TAES-SBIC", "<7>-HVUL", "<3>-<4>", "<3>-BDIS", "LPER-<9>", "SBIC-HVUL", "TURA-HVUL", "ATAU-BDIS", "SITA-SITA", "ZMAY-HVUL", "SITA-<9>", "<1>-<7>", "OSAT-HVUL", "<2>-OSAT", "<5>-HVUL", "HVUL-<9>", "BDIS-BDIS", "BDIS-<9>", "TURA-TURA", "ZMAY-SBIC", "ZMAY-<7>", "OSAT-OSAT", "<2>-<7>", "LPER-LPER", "SBIC-SITA", "ATAU-ATAU", "<3>-SITA", "<1>-ZMAY", "LPER-HVUL", "SITA-HVUL", "<3>-<5>", "BDIS-HVUL", "<3>-LPER", "<1>-SBIC", "OSAT-LPER", "<4>-<8>", "OSAT-<5>", "<8>-BDIS", "SBIC-SBIC", "TURA-BDIS", "<7>-BDIS", "OSAT-TURA", "<7>-ATAU", "SBIC-BDIS", "<3>-OSAT", "<8>-ATAU", "<7>-TURA", "<8>-TURA", "TURA-LPER", "<5>-TURA", "SBIC-ZMAY", "SBIC-<7>", "<6>-SITA", "<2>-ZMAY", "<4>-SITA", "ZMAY-BDIS", "<5>-ATAU", "OSAT-ATAU", "SBIC-ATAU", "<7>-SITA", "<7>-<8>", "HVUL-ATAU", "SBIC-TURA", "ATAU-LPER", "OSAT-BDIS", "ZMAY-ATAU", "<4>-LPER", "<5>-BDIS", "ZMAY-SITA", "<8>-LPER", "<2>-SBIC", "ZMAY-TURA", "<3>-ZMAY", "LPER-BDIS", "HVUL-BDIS", "LPER-ATAU", "<7>-<7>", "HVUL-TURA", "<7>-LPER", "SITA-BDIS", "BDIS-ATAU", "<3>-SBIC", "SITA-ATAU", "TURA-OSAT", "SBIC-LPER", "<8>-<5>", "LPER-TURA", "LPER-OSAT", "LPER-<5>", "SBIC-<8>", "TURA-SITA", "SITA-TURA", "ZMAY-LPER", "BDIS-TURA", "<4>-OSAT", "ATAU-<2>", "<3>-<3>", "SITA-LPER", "BDIS-LPER", "HVUL-LPER", "<4>-ZMAY", "<5>-SITA", "ATAU-SITA", "ZMAY-<8>", "ATAU-OSAT", "<8>-OSAT", "<3>-<7>", "<4>-<6>", "<4>-<5>", "<4>-SBIC", "OSAT-SITA", "<6>-ZMAY", "<7>-OSAT", "<6>-SBIC", "TURA-<5>", "HVUL-SITA", "TURA-TAES", "LPER-SITA", "SBIC-OSAT", "BDIS-SITA", "<5>-<5>", "ZMAY-OSAT", "<5>-<8>", "<4>-<7>", "SITA-OSAT", "ATAU-<5>", "BDIS-OSAT", "HVUL-OSAT", "<6>-<7>", "TURA-ZMAY", "<7>-<5>", "ATAU-<3>", "<5>-ZMAY", "SBIC-<5>", "ATAU-ZMAY", "TURA-SBIC", "OSAT-ZMAY", "SITA-ZMAY", "SITA-SBIC", "BDIS-<5>", "ZMAY-<5>", "<5>-SBIC", "HVUL-ZMAY", "ATAU-SBIC", "LPER-ZMAY", "BDIS-ZMAY", "OSAT-SBIC", "HVUL-SBIC", "LPER-SBIC", "BDIS-SBIC", "SITA-<7>", "SITA-<8>", "SITA-<5>", "HVUL-<5>", "<8>-<8>", "TURA-<8>", "TURA-<7>", "<8>-<4>", "<5>-<7>", "OSAT-<8>", "ATAU-<7>", "OSAT-<7>", "ATAU-<8>", "TURA-<3>", "ATAU-<4>", "LPER-<8>", "HVUL-<7>", "OSAT-<4>", "LPER-<7>", "BDIS-<7>", "BDIS-<8>", "TURA-<4>", "HVUL-<8>", "<7>-<3>", "<7>-<4>", "<8>-<3>", "<7>-<6>", "OSAT-<3>", "SBIC-<3>", "<4>-<4>", "LPER-<4>", "<5>-<3>", "SBIC-<4>", "ZMAY-<3>", "<5>-<6>", "<5>-<4>", "BDIS-<4>", "OSAT-<6>", "BDIS-<3>", "HVUL-<2>", "HVUL-<3>", "ZMAY-<4>", "LPER-<3>", "HVUL-<4>", "SITA-<3>", "TURA-<2>", "TURA-<1>", "SBIC-<6>", "ZMAY-<6>", "TURA-<6>", "SITA-<4>", "LPER-<6>", "<7>-<2>", "ATAU-<6>", "<8>-<2>", "OSAT-<2>", "<5>-<2>", "<6>-<6>", "BDIS-<6>", "SBIC-<2>", "SITA-<6>", "ZMAY-<2>", "LPER-<2>", "HVUL-<6>", "BDIS-<2>", "SITA-<2>", "<1>-<1>", "TAES-TAES", "ATAU-<1>", "ATAU-TAES", "HVUL-<1>", "<7>-<1>", "<8>-<1>", "OSAT-<1>", "<5>-<1>", "SBIC-<1>", "ZMAY-<1>", "LPER-<1>", "BDIS-<1>", "SITA-<1>", "<7>-TAES", "HVUL-TAES", "<8>-TAES", "OSAT-TAES", "<5>-TAES", "SBIC-TAES", "ZMAY-TAES", "LPER-TAES", "BDIS-TAES", "SITA-TAES"]}], {"autosize": false, "title": "GRAMPA Results: wheat_all_D_out.txt", "paper_bgcolor": "#fffae6", "plot_bgcolor": "#e1e1ea", "yaxis": {"titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "Score"}, "height": 500, "width": 1000, "titlefont": {"family": "Arial, sans-serif", "size": 30}, "xaxis": {"tickangle": 90, "titlefont": {"color": "#737373", "family": "Arial, sans-serif", "size": 20}, "title": "H1-H2 Node"}, "margin": {"b": 150, "r": 20, "pad": 0, "t": 70, "l": 70}}, {"linkText": "Export to plot.ly", "showLink": true})</script>';


outfilename = "../../results/wheat/" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, p1=p1, p2=p2, p3=p3, footer=footer));