import grampa_lib.recontree as RT
import grampa_lib.mul_tree as MT

def orthLabel(gene_trees, min_maps, min_tree, min_clade):

	orthfile = open(globs.orth_file_name, "w");
	lfile = open(globs.labeled_tree_file, "w");
	orthfile.write(MT.mulPrint(min_tree, min_clade));
	orthfile.write("\n------------------\n");

	for gene_num in gene_trees:
		gt, ginfo = gene_trees[gene_num];
		# Retrieve gene tree info.

		cur_results = min_maps[gene_num];
		# Getting the maps for the current gene tree.

		if len(cur_results) != 1:
			orthfile.write(str(gene_num+1) + "\t* " + str(len(cur_results)) + " maps tied for lowest score. Mapping orthologies to each set of maps!\n");
			lfile.write(str(gene_num+1) + "\t* " + str(len(cur_results)) + " maps tied for lowest score. Labeling trees for each set of maps!\n");

		for each_result in cur_results:
			maps, dups = each_result[3], each_result[4];

			spec_genes = {};
			for gene in ginfo:
				if ginfo[gene][2] == 'tip':
					spec = gene[gene.rfind("_")+1:];
					if spec in min_clade:
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

			lfile.write(str(gene_num+1) + "\t" + gt + "\n");

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

			orthfile.write(outline[:-1] + "\n");

	lfile.close();
	orthfile.close();



