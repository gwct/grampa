############################################################
# For GRAMPA site, 12.19
# This generates the file "performance.html"
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
				<h1>Program performance</h1>

				<h4>This page details the run time and memory usage of GRAMPA on several datasets. This might help you decide how to allocate resources
					for your dataset.</h4>

				<h2>Comparisons between real datasets</h2>

				<div class="table-container">
					<table class="stripe-table" id="spec_table">
						<thead><tr><th>Dataset</th><th># Taxa</th><th># Gene trees</th><th># MUL-trees considered</th></tr></thead>
						<tbody>
							<tr><td><a href="results/yeast/yeast.html">Yeast</a></td><td>27</td><td>4,623</td><td>13</td></tr>
							<tr><td><a href="results/wheat/wheat_all.html">Wheat</a></td><td>10</td><td>9,147</td><td>304</td></tr>
						</tbody>
					</table>
				</div>

				<h4>These datasets were run on the <a href="https://kb.iu.edu/d/bbhh" target="_blank">Mason computing cluster</a> at Indiana University.</h4>


				<div class="row">
					<div class="col-12-24">
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/real_cpu.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>									
					</div>

					<div class="col-12-24">
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/real_ram.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>				
					</div>
				</div>

				<h2>Comparisons between operating systems</h2>

				<div class="table-container">
					<table class="stripe-table" id="spec_table">
						<thead><tr><th>Dataset</th><th># Taxa</th><th># Gene trees</th><th># MUL-trees considered</th></tr></thead>
						<tbody>
							<tr><td>Simulation 3a</td><td>7</td><td>1,000</td><td>128</td></tr>
						</tbody>
					</table>
				</div>

				<div class="row"><div class="col-1" id="divider_row"></div></div>

				<div class="table-container">
					<table class="stripe-table" id="spec_table">
					<thead><tr><th>OS</th><th>Processor</th><th>Memory (GB)</th></tr></thead>
					<tbody>
						<tr><td>Mac OS X Yosemite (10.10.3)</td><td>Intel Core i7 @ 3.1 GHz</td><td>16</td></tr>
						<tr><td>Windows 7 Ultimate (64 bit)</td><td>Intel Core i5-750 @ 2.67 GHz</td><td>12</td></tr>
						<tr><td>Red Hat Enterprise Linux 6.x (login node on <a href="https://kb.iu.edu/d/bbhh" target="_blank">Mason</a>)</td><td>Two Intel Xeon E5-2600</td><td>24</td></tr>
					</tbody>
					</table>
				</div>

				<div class="row"><div class="col-1" id="divider_row"></div></div>

				<div class="row">
					<div class="col-12-24">
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/sim_cpu.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>									
					</div>

					<div class="col-12-24">
						<div class="row img-row">
							<div class="col-5-24 img-margin-left"></div>
							<div class="col-14-24 img-col">						
								<img class="grid-img" src="img/sim_ram.png">			
							</div>
							<div class="col-5-24 img-margin-right"></div>
						</div>				
					</div>
				</div>
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
pagefile = "performance.html";
print("Generating " + pagefile + "...");
title = "GRAMPA - Performance"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile, "", "results/yeast/", "results/wheat/");
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));