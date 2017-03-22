import sys, re
sys.path.append("C:\Users\Gregg\Desktop\Box Sync\mul_recon_project\grampa_sims_ohno\grampa-sim-dev\lib");
import recontree as RT

infilename = sys.argv[1];
labfilename = sys.argv[2];
hybrid_spec = sys.argv[3];

hybrid_spec = hybrid_spec.split(",");

tree_list = [];
num_trees = 0;
for tree in open(infilename):
	if tree[0] == "#":
		tree_list.append([]);
		continue;
	tree = re.sub("<[\d]+>", "", tree);
	tinfo, t = RT.treeParse(tree.strip());
	tree_list.append([tinfo,t]);
	num_trees += 1;

print num_trees;

num_same = 0;

for line in open(labfilename):
	line = line.strip().split("\t");
	gt_id = int(line[0])-1;
	tree = re.sub("<[\d]+>", "", line[1]);
	#print gt_id, "\t", tree67yy;
	tinfo, t = RT.treeParse(tree.strip());

	orig_tinfo = tree_list[gt_id][0];
	orig_tree = tree_list[gt_id][1];

	carrots = [];
	plusses = [];

	for node in tinfo:
		if tinfo[node][2] != 'tip' or node[node.rfind("_")+1:].replace("^","").replace("+","") not in hybrid_spec:
			continue;

		#print node;

		if "^" in node:
			carrots.append(node.replace("^", ""));
		if "+" in node:
			plusses.append(node.replace("+",""));



	stars = [];
	norms = [];

	for node in orig_tinfo:
		if orig_tinfo[node][2] != 'tip' or node[node.rfind("_")+1:].replace("*","") not in hybrid_spec:
			continue;

		#print node;

		if "*" in node:
			stars.append(node.replace("*",""));
		else:
			norms.append(node);



	if (set(norms) == set(carrots) and set(stars) == set(plusses)) or (set(norms) == set(plusses) and set(stars) == set(carrots)):
		num_same += 1;
	else:
		print t;
		print "CARROTS:", carrots;
		print "PLUSSES:", plusses;
		print orig_tree;
		print "STARS:", stars;
		print "NORMS:", norms;
		print "---------";
print num_same;