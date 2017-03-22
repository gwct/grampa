## Given a GT, a MT, an LCA in the MT, and a set of maps
import sys, os, re, time, recontree as RT, reconcore as RC, mul_recon as ALG

def lcaCheck(gt, ginfo, maps, dups, hybrid_clade, gene_num, loutfile, orthoutfile):

	spec_genes = {};
	for gene in ginfo:
		if ginfo[gene][2] == 'tip':
			spec = gene[gene.rfind("_")+1:];
			if spec in hybrid_clade:
				cur_map = maps[gene];
				if "*" in cur_map[0]:
					labeled_gene = gene + "+";
					gt = gt.replace(gene, labeled_gene);
				else:
					labeled_gene = gene + "^";
					gt = gt.replace(gene, labeled_gene);

				if spec not in spec_genes:
					spec_genes[spec] = [gene];
				else:
					spec_genes[spec].append(gene);

	loutfile.write(str(gene_num+1) + "\t" + gt + "\n");


	flag = 0;
	done = [];

	outline = str(gene_num+1) + "\t";

	for spec in spec_genes:
		if len(spec_genes[spec]) == 1:
			outline += spec_genes[spec][0] + "-SINGLE\t";

		else:
			for gene1 in spec_genes[spec]:
				for gene2 in spec_genes[spec]:
					if gene2 == gene1 or [gene1,gene2] in done or [gene2,gene1] in done:
						continue;
					done.append([gene1,gene2]);
					outline += gene1 + "-" + gene2 + "-";
					cur_lca = RT.LCA([gene1,gene2],ginfo)[0];
					if dups[cur_lca] == 1:
						outline += "PARALOG\t";
					else:
						outline += "HOMOEOLOG\t";

	orthoutfile.write(outline[:-1] + "\n");
	# orthoutfile.write(gt + "\n");
	# orthoutfile.write(str(maps) + "\n");
	# orthoutfile.write(str(dups) + "\n");
	# orthoutfile.write("---------------\n");

	# for spec in spec_genes:
	# 	specoutname = out_file_prefix + "_" + spec + ".txt";
	# 	specout = open(specoutname, "a");

	# 	#print spec, ":", spec_genes[spec];
	# 	if len(spec_genes[spec]) > 1:
	# 		flag = 1;
	# 	outline = str(gene_num+1) + "\t";
	# 	if len(spec_genes[spec]) == 1:
	# 		outline + spec_genes[spec][0] + "-SINGLE\t";

	# 	else:

	# 		for gene1 in spec_genes[spec]:
	# 			for gene2 in spec_genes[spec]:
	# 				if gene2 == gene1 or [gene1,gene2] in done or [gene2,gene1] in done:
	# 					continue;
	# 				done.append([gene1,gene2]);
					
	# 				outline += gene1 + "-" + gene2 + "-";

	# 				cur_lca = RT.LCA([gene1,gene2],ginfo)[0];
	# 				outline += gene1 + "-" + gene2 + "-" + cur_lca + "-";

	# 				if dups[cur_lca] == 1:
	# 					outline += "PARALOG\t";
	# 				else:
	# 					outline += "HOMOEOLOG\t";

	# 				# print outline;













	# spec_genes = {};
	# for gene in ginfo:
	# 	if ginfo[gene][2] == 'tip':
	# 		spec = gene[gene.rfind("_")+1:];
	# 		if spec in hybrid_clade:
	# 			cur_map = maps[gene];
	# 			if "*" in cur_map[0]:
	# 				labeled_gene = gene + "+";
	# 				gt = gt.replace(gene, labeled_gene);
	# 			else:
	# 				labeled_gene = gene + "^";
	# 				gt = gt.replace(gene, labeled_gene);

	# 			if spec not in spec_genes:
	# 				spec_genes[spec] = [labeled_gene];
	# 			else:
	# 				spec_genes[spec].append(labeled_gene);

	# loutfile.write(str(gene_num+1) + "\t" + gt + "\n");

	# flag = 0;
	# done = [];

	# outline = str(gene_num+1) + "\t";

	# for spec in spec_genes:
	# 	if len(spec_genes[spec]) == 1:
	# 		outline += spec_genes[spec][0] + "-SINGLE\t";

	# 	else:
	# 		for gene1 in spec_genes[spec]:
	# 			for gene2 in spec_genes[spec]:
	# 				if gene2 == gene1 or [gene1,gene2] in done or [gene2,gene1] in done:
	# 					continue;
	# 				done.append([gene1,gene2]);

	# 				# spec1 = gene1[gene1.index("_")+1:];
	# 				# spec2 = gene2[gene2.index("_")+1:];

	# 				outline += gene1 + "-" + gene2 + "-";

	# 				# if spec1 == spec2:
	# 				if len([g for g in [gene1,gene2] if "+" in g]) == 2 or len([g for g in [gene1,gene2] if "^" in g]) == 2:
	# 					outline += "PARALOG\t";
	# 				else:
	# 					outline += "HOMOEOLOG\t";

	# 				# else:
	# 				# 	if len([g for g in [gene1,gene2] if "+" in g]) == 2 or len([g for g in [gene1,gene2] if "^" in g]) == 2:
	# 				# 		outline += "PARALOG\t";
	# 				# 	else:
	# 				# 		outline += "ORTHOLOG\t";				

	# 				# print outline;
	# orthoutfile.write(outline[:-1] + "\n");

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

	#print "GT:", gt;
	#print "MT:", mt;

	#print "LCA DUPS:", dups;

	# log_types = {};
	# spec_genes = {};
	# for spec in hybrid_spec:
	# 	log_types[spec] = [];
	# 	spec_genes[spec] = [];

	# for gene in maps:
	# 	if ginfo[gene][2] == 'tip':
	# 		spec = gene[gene.index("_")+1:];
	# 		if spec in hybrid_spec:
	# 			spec_genes[spec].append(gene);

	# flag = 0;
	# done = [];

	# for spec in spec_genes:
	# 	specoutname = out_file_prefix + "_" + spec + ".txt";
	# 	specout = open(specoutname, "a");

	# 	#print spec, ":", spec_genes[spec];
	# 	if len(spec_genes[spec]) > 1:
	# 		flag = 1;
	# 	outline = str(gene_num+1) + "\t";
	# 	if len(spec_genes[spec]) == 1:
	# 		outline + spec_genes[spec][0] + "-SINGLE\t";

	# 	else:

	# 		for gene1 in spec_genes[spec]:
	# 			for gene2 in spec_genes[spec]:
	# 				if gene2 == gene1 or [gene1,gene2] in done or [gene2,gene1] in done:
	# 					continue;

	# 				done.append([gene1,gene2]);

	# 				cur_lca = RT.LCA([gene1,gene2],ginfo)[0];
	# 				outline += gene1 + "-" + gene2 + "-" + cur_lca + "-";

	# 				if dups[cur_lca] == 1:
	# 					outline += "PARALOG\t";
	# 				else:
	# 					outline += "HOMOEOLOG\t";

	# 				# print outline;
	# 	specout.write(outline[:-1] + "\n");

	# 	specout.close();

	# if (gene_num+1) == 193:
	# 	print "GT:", gt;
	# 	print "MT:", mt;
	# 	print "MAPS:", maps;
	# 	print "DUPS:", dups;
	# 	sys.exit();


	#print "=========================================";


	#if gene_num > 20:
	#	sys.exit();






	# singles = [];
	# for gene in maps:
	# 	flag = 0;
	# 	if ginfo[gene][2] == 'tip':
	# 		spec = gene[gene.index("_")+1:];
	# 		if spec in hybrid_spec:
	# 			for group in gt_groups:
	# 				if gene in group:
	# 					flag = 1;
	# 					break;
	# 			if flag == 0:
	# 				singles.append(gene);

	# for each in singles:
	# 	gt_groups.append([each]);

	# star_map = 0;
	# norm_map = 0;

	# for group in gt_groups:
	# 	if "*" in maps[group[0]][0]:
	# 		star_map += 1;
	# 	else:
	# 		norm_map += 1;

	# return star_map, norm_map, len(gt_groups);




	# found = 0;
	# star_spec_count = 0;
	# norm_spec_count = 0;
	# both_spec_count = 0;
	# spec_miss_count = 0;

	# for spec in hybrid_spec:
	# 	gene_copies = [node for node in ginfo if ginfo[node][2] == 'tip' and "_" + spec in node];
	# 	#print gene_copies;

	# 	if gene_copies == []:
	# 		spec_miss_count += 1;
	# 		continue;

	# 	gene_copy_maps = [maps[m][0] for m in maps if m in gene_copies];

	# 	num_star = len([m for m in gene_copy_maps if "*" in m]);
	# 	if num_star == len(gene_copy_maps):
	# 		star_spec_count += 1;
	# 	elif num_star == 0:
	# 		norm_spec_count += 1;
	# 	else:
	# 		both_spec_count += 1;

	# return star_spec_count, norm_spec_count, both_spec_count, spec_miss_count;
	# COUNT THE NUMBER OF SPECIES MAPPING TO STAR, NORMAL, BOTH, OR MISSING



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