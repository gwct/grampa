############################################################
# For GRAMPA site, 12.19
# This generates the file "index.html"
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

				<div class="row img-row img-no-mobile">
					<div class="col-8-24 img-margin-left"></div>
					<div class="col-8-24 img-col">
						<img class="grid-img" id="logo_main" src="img/logo2.png">
					</div>
					<div class="col-8-24 img-margin-right"></div>
				</div>

				<h1>
					GRAMPA: Gene-tree Reconciliation Algorithm with MUL-trees for Polyploid Analysis
				</h1>

				<div class="row" id="title-row">
					<div class="col-2-24 margin"></div>
					<div class="col-20-24" id="title-col">
						<center>
							<a href="https://bioconda.github.io/recipes/grampa/README.html" target="_blank">
								<img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" />
								<img src="https://anaconda.org/bioconda/grampa/badges/platforms.svg" />
								<img src="https://img.shields.io/conda/vn/bioconda/grampa?label=version" />
								<img src="https://anaconda.org/bioconda/grampa/badges/latest_release_date.svg" />
								<img src="https://img.shields.io/conda/dn/bioconda/grampa.svg?style=flat" />
							</a>
							<a href="https://github.com/gwct/grampa" target="_blank">
								<img src="https://img.shields.io/github/commits-since/gwct/grampa/v1.4.0" />
								<img alt="GitHub" src="https://img.shields.io/github/license/gwct/grampa">
							</a>
							
						</center>
					</div>
					<div class="col-2-24 margin"></div>
				</div>

				<p>
					GRAMPA uses gene tree topologies for inferring the presence and mode of polyploidy in a set of species, placing a whole genome duplication event on 
					a phylogeny, and counting gene duplications and losses in the presence of polyploidy.
				</p>

				<h2>Install</h2>

				<p>
					<a href="https://bioconda.github.io/recipes/grampa/README.html" target="_blank">Install with bioconda</a>:
				</p>

				<div class="row">
					<div class="col-7-24 margin"></div>
					<div class="col-10-24" id="install">
						<center><pre class="cmd"><code>conda install -c bioconda grampa</code></pre></center>
					</div>
					<div class="col-7-24 margin"></div>
				</div>	

				<p>
					Otherwise, just download the program straight <a href="https://github.com/gwct/grampa/releases/latest" target="_blank">from github</a> by clicking 
					the download button below. The program is just a python script with no other dependencies, so it should work out of the box!
				</p>

				<h2>Citation</h2>
				<p id="paper">
					Thomas GWC, Ather SH, Hahn MW. 2017. Gene tree reconciliation with MUL-trees to resolve polyploid analysis. <em>Systematic Biology.</em>
					66(6):1007-1018. <a href="https://doi.org/10.1093/sysbio/syx044" target="_blank">10.1093/sysbio/syx044</a>
				</p>
				<div class="row" id="buttons_container">
					<div class="col-2-24 button">
						<a class="button-secondary id="install_btn" href="about.html">About&nbsp;&raquo;</a>
					</div>

					<div class="col-1-24 btn-sep"></div>

					<div class="col-2-24 button">
						<a class="button-secondary id="install_btn" href="readme.html">README&nbsp;&raquo;</a>
					</div>

					<div class="col-1-24 btn-sep"></div>

					<div class="col-2-24 button">
						<a class="button-secondary id="install_btn" href="https://github.com/gwct/grampa/releases/latest" target="_blank">Download&nbsp;&raquo;</a>
					</div>

					<div class="col-16-24 btn-margin"></div>
				</div>
			</div>
		</div>
		<div class="col-3-24" id="margin"></div>
	</div>

	<div class="row" id="task_row">
		<div class="col-3-24" id="margin"></div>
		<div class="col-18-24" id="task_col">
			<div class="row" id="task_row_inner">
				<div class="col-1-3" id="task">
					<h2>Identify the mode of polyploidy</h2>
					<p id="ex_p">GRAMPA may be able to differentiate between the cases of allopolyploidy, autopolyploidy, and no polyploidy.</p>
				</div>
				<div class="col-1-3" id="task">
					<h2>Place a WGD on a phylogeny</h2>
					<p id="ex_p">If you know a WGD has taken place on your phylogeny, GRAMPA can place it on the tree.</p>
				</div>
				<div class="col-1-3" id="task">
					<h2>Count duplications and losses</h2>
					<p id="ex_p">GRAMPA can accurately count duplications and losses, regardless of the presence of polyploidy.</p>
				</div>
			</div>
			<div class="row" id="btn_row_inner">
				<div class="col-2-24 button">
					<a class="button-secondary" href="example1.html">Example&nbsp&raquo;</a>
				</div>
				<div class="col-6-24"></div>

				<div class="col-2-24 button">
					<a class="button-secondary" href="example2.html">Example&nbsp&raquo;</a>
				</div>
				<div class="col-6-24"></div>

				<div class="col-2-24 button">
					<a class="button-secondary" href="example3.html">Example&nbsp&raquo;</a>
				</div>								
				<div class="col-6-24"></div>
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
pagefile = "index.html";
print("Generating " + pagefile + "...");
title = "GRAMPA"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));