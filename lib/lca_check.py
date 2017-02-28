## Given a GT, a MT, an LCA in the MT, and a set of maps
import sys, os, re, time, recontree as RT, reconcore as RC, mul_recon as ALG

def lcaCheck(mt, minfo, gt, ginfo, maps, lca, hybrid_spec):
	# mt = "((A,(B,C)),((D,B*),E))"
	# gt = "((1_A,(1_B,1_C)),(2_B,1_E))"
	# gt = "(((1_A,1_B),(2_B,1_C)),(1_D,1_E))"

	# minfo, mt = RT.treeParse(mt);
	# ginfo, gt = RT.treeParse(gt);

	# print mt;	
	# print gt;

	# maps = {
	# 	'1_A' : ['A'],
	# 	'1_B' : ['B'],
	# 	'1_C' : ['C'],
	# 	'2_B' : ['B'],
	# 	'1_E' : ['E'],
	# 	'<1>' : ['<1>'],
	# 	'<2>' : ['<2>'],
	# 	'<3>' : ['<4>'],
	# 	'<4>' : ['<5>']
	# }

	# lca = "<5>";

	# hybrid_spec = ['B'];

	found = 0;
	star_spec_count = 0;
	norm_spec_count = 0;
	both_spec_count = 0;
	spec_miss_count = 0;

	for spec in hybrid_spec:
		gene_copies = [node for node in ginfo if ginfo[node][2] == 'tip' and "_" + spec in node];
		#print gene_copies;

		if gene_copies == []:
			spec_miss_count += 1;
			continue;

		gene_copy_maps = [maps[m][0] for m in maps if m in gene_copies];

		num_star = len([m for m in gene_copy_maps if "*" in m]);
		if num_star == len(gene_copy_maps):
			star_spec_count += 1;
		elif num_star == 0:
			norm_spec_count += 1;
		else:
			both_spec_count += 1;

	return star_spec_count, norm_spec_count, both_spec_count, spec_miss_count;




		# done = [];

		# for copy1 in gene_copies:
		# 	for copy2 in gene_copies:
		# 		if copy2 == copy1 or [copy1,copy2] in done or [copy2,copy1] in done:
		# 			continue;

		# 		done.append([copy1,copy2]);

		# 		#print copy1;
		# 		#print copy2;

		# 		gene_lca = RT.LCA([copy1,copy2],ginfo)[0];
		# 		#print gene_lca;

		# 		lca_map = maps[gene_lca][0];
		# 		#print lca_map;

		# 		if lca_map == lca:
		# 			found = 1;
		# CHECKS FOR GENES WITH LCA MAPPING TO LCA IN MUL TREE

	# return found;