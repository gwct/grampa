import reconcore as RC, mul_recon as ALG, gene_tree as GT

def stdLCA(st, sinfo, gene_trees_filtered, outfilename, detoutfilename, lca_opt, maps_opt, v):
# Only do standard recon if the input tree is a singly-labeled tree.
	st_score = 0;
	gene_num = -1;
	gt_maps = [];
	gt_dups = [];

	tot_node_counts = {};
	for node in sinfo:
		tot_node_counts[node] = [0,0];

	RC.printWrite(detoutfilename, v, "# ---------------------------");
	RC.printWrite(detoutfilename, v, "ST:" + st);
	for gene_tree in gene_trees_filtered:
		gene_num += 1;
		if len(gene_tree) == 1:
			gt_maps.append('');
			gt_dups.append('');
			continue;
		# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

		outline = "GT-" + str(gene_num+1) + " to ST\t";
		gt, ginfo = gene_tree;
		# Retrieves the gene tree info.

		maps = {};
		for g in ginfo:
			if ginfo[g][2] == 'tip':
				speclabel = g[g.rfind("_")+1:];
				maps[g] = [speclabel];
			else:
				maps[g] = [];
		# Initialize the maps.

		st_maps, st_dups, st_loss = ALG.reconLCA(ginfo, sinfo, maps);
		gt_maps.append(st_maps);
		gt_dups.append(st_dups);

		num_dups = sum(st_dups.values());
		num_loss = sum(st_loss.values());

		st_mut_score = num_dups + num_loss;
		st_score += st_mut_score;
		outline = outline + str(num_dups) + "\t" + str(num_loss) + "\t" + str(st_mut_score);
		RC.printWrite(detoutfilename, v, outline);
		if maps_opt:
			GT.detailedOut(st_maps, st_dups, st_loss, v, detoutfilename);
		# Call the recon algorithm and aggregate scores.

	RC.printWrite(detoutfilename, v, "Total parsimony score for ST: " + str(st_score));

	RC.printWrite(detoutfilename, v, "# ---------------------------");
	RC.printWrite(outfilename, 0, "ST\t\t\t" + st + "\t" + str(st_score));

 	# Print the total score and total branch scores for the singly-labeled tree.
	### End standard recon block.
	lowest_score = ["ST", st, sinfo, "", gt_maps, gt_dups, st_score];
	return lowest_score;