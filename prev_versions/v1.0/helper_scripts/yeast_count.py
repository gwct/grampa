import sys, bar_plot

infilename = "../../yeast_test/yeast_test_orthologies.txt";
hybrid_spec = ['51660A','51914A','588726A','1071379A','CANGA','KAZAF','NAUDC','SACBA','SACCA','TETPH','VANPO','YEAST'];

single_counts = {};
para_counts = {};
homoeo_counts = {};

spec_single_counts = { spec : 0 for spec in hybrid_spec };
spec_orth_counts = {};
done = [];
for spec1 in hybrid_spec:
	for spec2 in hybrid_spec:
		if spec1 + "-" + spec2 in spec_orth_counts or spec2 + "-" + spec1 in spec_orth_counts:
			continue;
		spec_orth_counts[spec1 + "-" + spec2] = [0,0,0];

#print spec_orth_counts;

for line in open(infilename):
	line = line.strip().split("\t");
	gt_id = line[0];
	genes = line[1:];

	gt_singles = 0;
	gt_paras = 0;
	gt_homoeos = 0;

	for gene in genes:
		gene = gene.replace("^","").replace("+","");
		gene = gene.split("-");

		#print gene;
		if len(gene) == 2:
			spec = gene[0][gene[0].index("_")+1:];
			# print spec;

			spec_single_counts[spec] += 1;
			gt_singles += 1;

		else:
			spec1 = gene[0][gene[0].index("_")+1:];
			spec2 = gene[1][gene[1].index("_")+1:];
			orth_type = gene[2];

			# print spec1;
			# print spec2;
			# print orth_type;

			if spec1 + "-" + spec2 in spec_orth_counts:
				key = spec1 + "-" + spec2;
			else:
				key = spec2 + "-" + spec1;

			if orth_type == "PARALOG":
				gt_paras += 1;
				spec_orth_counts[key][0] += 1;

			elif orth_type == "HOMOEOLOG":
				gt_homoeos += 1;
				spec_orth_counts[key][1] += 1;
	# print gt_singles;
	# print gt_paras;
	# print gt_homoeos;

	if gt_singles not in single_counts:
		single_counts[gt_singles] = 1;
	else:
		single_counts[gt_singles] += 1;

	if gt_paras not in para_counts:
		para_counts[gt_paras] = 1;
	else:
		para_counts[gt_paras] += 1;

	if gt_homoeos not in homoeo_counts:
		homoeo_counts[gt_homoeos] = 1;
	else:
		homoeo_counts[gt_homoeos] += 1;


p_genes = [];
p_counts = [];

for x in range(1,51):
	p_genes.append(x);
	if x in para_counts:
		p_counts.append(para_counts[x]);
	else:
		p_counts.append(0);

p_xtitle = "# of paralogs";
p_ytitle = "# of genes";
maintitle = "";
outname = "yeast_paralogs.html";

bar_plot.barPlot(p_genes,p_counts,p_xtitle,p_ytitle,maintitle,outname,barcol='rgb(0,102,51)',plotcol='#e1e1ea',bgcol='#fffae6',w=1000,h=500,bmar=150);

h_genes = [];
h_counts = [];

for x in range(1,51):
	h_genes.append(x);
	if x in homoeo_counts:
		h_counts.append(homoeo_counts[x]);
	else:
		h_counts.append(0);

h_xtitle = "# of homoeologs";
h_ytitle = "# of genes";
maintitle = "";
outname = "yeast_homoeologs.html";

bar_plot.barPlot(h_genes,h_counts,h_xtitle,h_ytitle,maintitle,outname,barcol='rgb(0,102,51)',plotcol='#e1e1ea',bgcol='#fffae6',w=1000,h=500,bmar=150);

outname = "side.html";
data = [[p_genes,p_counts,"Paralogs"],[h_genes,h_counts,"Homoeologs"]];
bar_plot.sideBarPlot(data,maintitle,outname,barcol='rgb(0,102,51)',plotcol='#e1e1ea',bgcol='#fffae6',w=1000,h=500,bmar=150);

# for c in homoeo_counts:
# 	print c, "\t", homoeo_counts[c];








