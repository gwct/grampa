#!/usr/bin/python
#############################################################################
# The main algorithmic functions for MUL-reconciliation mapping.
# Gregg Thomas
# Fall 2015, Combo algorithm implemented Spring 2016
#############################################################################

import sys, re, recontree as RT, mul_tree as MT, itertools

#############################################################################

def reconLCA(lca_ginfo, sinfo, lca_maps):
# The LCA reconciliation mapping algorithm.

	dups = {};
	losses = {};
	internal_nodes = [];

	internal_nodes = RT.sortNodes(lca_ginfo);
	# Sort the internal nodes for a post order traversal.

	for g in lca_ginfo:
		dups[g] = 0;
		losses[g] = 0;

	for g in internal_nodes:
		g = "<" + str(g) + ">";
		d1, d2 = RT.getDesc(g, lca_ginfo);
		is_dup = 0;
		g_clade = RT.getClade(g, lca_ginfo);
		clade_maps = [];
		for g_tip in g_clade:
			clade_maps.append(lca_maps[g_tip][0]);
		# Get the species in the clade of the current node. Then get all
		# the possible maps from those species.

		lca_maps[g].append(RT.LCA(clade_maps,sinfo)[0]);

		if lca_maps[g][0] == lca_maps[d1][0] or lca_maps[g][0] == lca_maps[d2][0]:
			dups[g] += 1;
			is_dup = 1;
		#Now, if the map of g is identical to one of its descendants, it is a duplication node.

		cur_depth = len(RT.nodeDepth(lca_maps[g][0],sinfo))

		if lca_ginfo[g][2] == 'root':
			losses[g] += cur_depth;
		# The number of losses at the root of the gene tree is equal to the depth of its map.

		d1_depth = len(RT.nodeDepth(lca_maps[d1][0],sinfo))
		d1_loss = (d1_depth - cur_depth - 1) + is_dup;
		losses[d1] += d1_loss;
		d2_depth = len(RT.nodeDepth(lca_maps[d2][0],sinfo))
		d2_loss = (d2_depth - cur_depth - 1) + is_dup;
		losses[d2] += d2_loss;
		# Counting losses for each of the descendents of the current node.

	return lca_maps, dups, losses;
	# Return the total number of duplication nodes.

#############################################################################
def getSis(gs_node, check_node, check_clade, gs_dict):
# Gets the hybrid and copy sister species.

	d1, d2 = RT.getDesc(gs_node, gs_dict);
	if d1 == check_node:
		sis_node = d2;
	elif d2 == check_node:
		sis_node = d1;

	sis_clade = RT.getClade(sis_node, gs_dict);
	if any(c in check_clade for c in sis_clade):
		return [];
	else:
		return sis_clade;

#############################################################################
def collapseGroups(mul_dict, gene_trees_filtered, checkfilename, spec_type, num_skipped, cap, v):
# The collapseGroups function goes through all gene tree-MUL-tree combos to collapse the groups. It also filters out
# gene trees that have more groups than the cap in any combination.
	checkfile = open(checkfilename, "a");
	for mul_num in mul_dict:
		gt_groups = [[] for gt in gene_trees_filtered];

		mt, minfo, hybrid_clade, hybrid_node, copy_node = mul_dict[mul_num][0], mul_dict[mul_num][1], mul_dict[mul_num][2], mul_dict[mul_num][3], mul_dict[mul_num][4]

		if spec_type == 's':
			checkfile.write("MT-" + str(mul_num) + ":" + MT.mulPrint(mt, hybrid_clade) + "\tH1 Node:" + hybrid_node + "\tH2 Node:" + copy_node + "\n");
		elif spec_type == 'm':
			checkfile.write("MT-" + str(mul_num) + ":" + MT.mulPrint(mt, hybrid_clade) + "\n");

		gene_num = -1;
		for gene_tree in gene_trees_filtered:
			gene_num += 1;
			if len(gene_tree) == 1:
				continue;
			# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

			gt,ginfo = gene_tree;
			outline = "GT-" + str(gene_num+1) + " to MT-" + str(mul_num) + "\t";
			
			internal_nodes = RT.sortNodes(ginfo);
			# Sort the internal nodes for a post order traversal.

			singles, groups = {}, {};

			for g in ginfo:
				if ginfo[g][2] == 'tip':
					if g[g.rfind("_")+1:] in hybrid_clade:
						cur_anc = ginfo[g][1];
						anc_clade = RT.getClade(cur_anc, ginfo);
						anc_clade.remove(g);
						singles[g] = anc_clade;
			# First, all hybrid species nodes in the gene tree are added to the singles list.
			## GETS SINGLETONS

			for g in internal_nodes:
				g = "<" + str(g) + ">";
				# Next, for any non-tip node, we find out if the species that define it can be grouped

				d1, d2 = RT.getDesc(g, ginfo);
				d1_clade = RT.getClade(d1, ginfo);
				d1_spec_clade = [spec[spec.rfind("_")+1:] for spec in d1_clade];
				d2_clade = RT.getClade(d2,ginfo);
				d2_spec_clade = [spec[spec.rfind("_")+1:] for spec in d2_clade];
				# The clades for the descendants of both nodes are retrieved, and their corresponding
				# species are stored.

				if all(s in hybrid_clade for s in d1_spec_clade) and all(s in hybrid_clade for s in d2_spec_clade):
				# If the descendants from both nodes are all hybrid clade species, then we may be able to group them.
					if not any(s in d2_spec_clade for s in d1_spec_clade):
					# However, only if there is not more than one copy of a species among the clades can they be grouped.
						cur_clade = RT.getClade(g, ginfo);
						cur_anc = ginfo[g][1];
						anc_clade = RT.getClade(cur_anc, ginfo);
						anc_clade = [spec for spec in anc_clade if spec not in cur_clade];

						cur_nodes = RT.getCladeNode(g, ginfo);
						for node in cur_nodes:
							if node in groups:
								del groups[node];

						groups[g] = [cur_clade, anc_clade];
				## CHECKS GROUPINGS

			for group in groups:
				for g in groups[group][0]:
					if g in singles:
						del singles[g];
			# Removes any singles that are in a group.

			final_groups = [];
			for node in groups:
				final_groups.append(groups[node]);
			for single in singles:
				final_groups.append([[single], singles[single]]);
			# Restructures the final groups and adds singles.

			sisters = {};

			if spec_type == 's':
				mul_hybrid_node = [n for n in minfo if set(RT.getClade(n, minfo)) == set(hybrid_clade)][0];
				copy_clade = [c + "*" for c in hybrid_clade];
				mul_copy_node = [n for n in minfo if set(RT.getClade(n, minfo)) == set(copy_clade)][0];
				# The copy clade is defined.
			elif spec_type == 'm':
				copy_clade = RT.getClade(copy_node, minfo);
				mul_hybrid_node = hybrid_node;
				mul_copy_node = copy_node;

			hybrid_anc = minfo[mul_hybrid_node][1];
			copy_anc = minfo[mul_copy_node][1];

			sisters[''] = getSis(hybrid_anc, mul_hybrid_node, copy_clade, minfo);
			sisters['*'] = getSis(copy_anc, mul_copy_node, hybrid_clade, minfo);
			# These lines get any sister species from the hybrid and copy clades in the MUL-tree and that
			# clade's corresponding map. If there are no sisters, it stores an empty list.

			groups, fixed_groups = [], [];

			for group in final_groups:
				group_sis = [spec[spec.rfind("_")+1:] for spec in group[1]];
				if group_sis == []:
					groups.append(group[0]);
					continue;

				if all(spec in sisters[''] for spec in group_sis):
					fixed_groups.append([group[0],'']);
				elif all(spec in sisters['*'] for spec in group_sis):
					fixed_groups.append([group[0],'*']);
				else:
					groups.append(group[0]);
			# This checks the sister species of all the groups for the gene tree. If all the sister species
			# of a group are also in the sister species of the hybrid or copy clade in the MUL-tree, then we
			# can fix the mapping of that node.
			## FINDS FIXED SISTER GROUPS

			gt_groups[gene_num] = [groups, fixed_groups];
			# Adding the groups and fixed groups to the current gt_groups.

			#num_singles = 
			num_groups = len(groups);
			num_fixed = len(fixed_groups);
			outline += str(num_groups) + "\t" + str(num_fixed) + "\t" + str(2**num_groups);
			if num_groups > cap:
				gene_trees_filtered[gene_num] = ["# Number of groups over group cap (-p set to " + str(cap) + ") -- Filtering."];
			 	outline += "\tNumber of groups over group cap (-p set to " + str(cap) + ") -- Filtering.";
			 	num_skipped += 1;
			checkfile.write(outline + "\n");
			# The call of the reconciliation algorithm! On the current gene tree with the current MUL-tree.

		mul_dict[mul_num][6] = gt_groups;
		checkfile.write("# ----------------------------------\n");
	checkfile.close();

	return mul_dict, gene_trees_filtered, num_skipped;

#############################################################################

def mulRecon(hybrid_clade, mt, minfo, gt, ginfo, cur_groups, cur_fixed, v):
# The basis of the MUL-reconciliation algorithm is that there are now nodes that
# have more than one possible map. We try all combinations of mappings for these
# nodes and find which combination results in the most parsimonious mutation score
# (# duplication + # losses).
#
# A few prelminary steps are taken to ensure the quickest mapping groups:
# 	1.  Identify whether the hybrid or copy clade in the MUL-tree have sister groups. If so, we can use
#		them to fix similar nodes in the gene tree.
#	2.  Find nodes that contain only one or zero copies of the hybrid node species and species from one
#		of the sister groups. Fix the mappings of these nodes.
#	3.  Any other nodes that contain only one or zero copies of the hybrid node species can be grouped
#		and should be mapped consistently, though we will still have to try both maps.
#	4.  Remaining single hybrid nodes must be tried with both maps.
#
# Once these steps are done, a list of node groups is obtained, for which we generate every combination
# of map and try to reconcile to the MUL-tree. A score is obtained for each combination and the minimum
# score is kept as the correct map.

	min_score = 99999;
	# This is the variable we are minimizing: score is defined as # duplications + # losses for any given
	# reconciliation.

	all_maps = [];
	lowest_maps = [];

	num_groups = len(cur_groups);

	#combo_ind = list(itertools.product(['','*'], repeat=len(node_ind)));
	#if v == -2:
	#	print "num combos", len(combo_ind);
	# We get all combinations of mappings for each node group. This is the time constraining step.
	map_num = 1;

	#combos = list(itertools.product(['','*'], repeat=len(node_ind)));
	for combo in itertools.product(['','*'], repeat=num_groups):
		group_map = [];
		for i in range(len(combo)):
			for node in cur_groups[i]:
				group_map.append(node + combo[i]);
		# This builds the current map for each group.

		for fixed in cur_fixed:
			for node in fixed[0]:
				group_map.append(node + fixed[1]);
		# This adds the fixed maps onto the current combination of group mappings.

		# Now we do LCA mapping for the current combination of maps for the hybrid clade species.
		maps = {};
		for g in ginfo:
			if ginfo[g][2] == 'tip':
				speclabel = g[g.rfind("_")+1:];
				if g in group_map:
					maps[g] = [speclabel];
				# If the node is in a hybrid clade, use the map in the current combination.
				elif g + "*" in group_map:
					maps[g] = [speclabel + "*"];
				else:
					maps[g] = [speclabel];
				# Otherwise, the map is just the species label.

			else:
				maps[g] = [];
			# And if the node is not a tip, the map is empty.
		# Initialization of the mapping dictionary, using the current combination of hybrid clade maps. 
		# [key]:[value] -> [gene tree node]:[map]

		maps, node_dups, node_loss = reconLCA(ginfo, minfo, maps);
		# Once the current maps have been initialized, we can simply call the normal LCA mapping algorithm

		num_dups = sum(node_dups.values());
		num_loss = sum(node_loss.values());

		cur_score = num_dups + num_loss;
 
		if cur_score <= min_score:
			min_score = cur_score;
		all_maps.append([cur_score, num_dups, num_loss, maps, node_dups, node_loss]);
			#min_node_counts = mul_node_counts;
		# If the current total score is smaller than the previous min, save it as the current min.

	for cur_map in all_maps:
		if cur_map[0] == min_score:
			lowest_maps.append(cur_map);

	return lowest_maps;
	# Return the maps with the minimal mutation score (dup + loss score).

	







