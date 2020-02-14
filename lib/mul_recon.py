#!/usr/bin/python
#############################################################################
# The main algorithmic functions for MUL-reconciliation mapping.
# Gregg Thomas
# Fall 2015, Combo algorithm implemented Spring 2016
#############################################################################

import os, itertools, recontree as RT, mul_tree as MT, reconcore as RC, gene_tree as GT, global_vars as globs
import pickle

#############################################################################

def reconLCA(lca_ginfo, sinfo, lca_maps, retmap=False):
# The LCA reconciliation mapping algorithm.

	internal_nodes = RT.sortNodes(lca_ginfo);
	# Sort the internal nodes for a post order traversal.

	score = 0;

	if retmap:
		dups, losses = {}, {};
		for g in lca_ginfo:
			dups[g], losses[g] = 0, 0;

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
			if retmap:
				dups[g] += 1;
			score += 1;
			is_dup = 1;
		#Now, if the map of g is identical to one of its descendants, it is a duplication node.

		cur_depth = len(RT.nodeDepth(lca_maps[g][0],sinfo))

		if lca_ginfo[g][2] == 'root':
			if retmap:
				losses[g] += cur_depth;
			score += cur_depth;
		# The number of losses at the root of the gene tree is equal to the depth of its map.

		d1_depth = len(RT.nodeDepth(lca_maps[d1][0],sinfo));
		d1_loss = (d1_depth - cur_depth - 1) + is_dup;
		score += d1_loss
		if retmap:
			losses[d1] += d1_loss;

		d2_depth = len(RT.nodeDepth(lca_maps[d2][0],sinfo))
		d2_loss = (d2_depth - cur_depth - 1) + is_dup;
		score += d2_loss;
		if retmap:
			losses[d2] += d2_loss;
		# Counting losses for each of the descendents of the current node.

	if retmap:
		return lca_maps, dups, losses;
	return score;
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
def collapseGroups(mul_input, gene_trees_filtered_cg, spec_type_cg, v, pickle_dir, nmt):
# The collapseGroups function goes through all gene tree-MUL-tree combos to collapse the groups.

	mul_num, mul_tree = mul_input;
	
	if v == 1:
		print("# " + RC.getDateTime() + " --> Collapsing groups for MUL-tree # " + str(mul_num) + " / " + str(nmt));

	if mul_num == 0:
		return mul_num, [];

	gt_groups = {};

	mt, minfo, hybrid_clade, hybrid_node, copy_node = mul_tree[0], mul_tree[1], mul_tree[2], mul_tree[3], mul_tree[4];

	for gene_num in gene_trees_filtered_cg:
		gene_tree = gene_trees_filtered_cg[gene_num];
		if len(gene_tree) == 1:
			continue;
		# If the gene tree was previously filtered, the list will only contain the filter message and it should be skipped here.

		gt,ginfo = gene_tree;
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

		if spec_type_cg == 's':
			mul_hybrid_node = [n for n in minfo if set(RT.getClade(n, minfo)) == set(hybrid_clade)][0];
			copy_clade = [c + "*" for c in hybrid_clade];
			mul_copy_node = [n for n in minfo if set(RT.getClade(n, minfo)) == set(copy_clade)][0];
			# The copy clade is defined.
		elif spec_type_cg == 'm':
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

	groupoutfile = os.path.join(pickle_dir, str(mul_num) + "_groups.pickle");
	pickle.dump(gt_groups, open(groupoutfile, "wb"));
	del groups, fixed_groups, final_groups, gene_trees_filtered_cg, gt_groups;

#############################################################################
def mulRecon(mul_input, gene_trees, v, pickle_dir, nmt, retmap=False):
# The basis of the MUL-reconciliation algorithm is that there are now nodes that
# have more than one possible map. We try all combinations of mappings for these
# nodes and find which combination(s) results in the most parsimonious mutation score
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
# Once these steps are done (in the collapseGroups function), a list of node groups is obtained, for
# which we generate every combination of map and try to reconcile to the MUL-tree. A score is obtained
# for each combination and the minimum score is kept as the correct map.

	mul_num, mul_tree = mul_input
	#main_output, det_output, min_num, min_score, min_maps, multiple_maps = {}, [], '', 9999999, {}, 0;
	# mulpicklefile = os.path.join(pickle_dir, str(mul_num) + "_tree.pickle");
	# mul_tree = pickle.load(open(mulpicklefile, "rb"));

	if v == 1:
		print("# " + RC.getDateTime() + " --> Reconciling to MUL-tree # " + str(mul_num) + " / " + str(nmt));

	min_maps = {};
	total_score = 0;

	if mul_num != 0:
		groupfilename = os.path.join(pickle_dir, str(mul_num) + "_groups.pickle");
		cur_groups = pickle.load(open(groupfilename, "rb"));

	for gene_num, gene_tree in gene_trees.items():
		gt, ginfo = gene_tree;

		gene_score = 99999;
		min_maps[gene_num] = [];

		if mul_num == 0:
			sinfo = mul_tree[1];

			init_maps = {};
			for g in ginfo:
				if ginfo[g][2] == 'tip':
					speclabel = g[g.rfind("_")+1:];
					init_maps[g] = [speclabel];
				else:
					init_maps[g] = [];
			# Initialize the maps.

			if retmap:
				maps, node_dups, node_loss = reconLCA(ginfo, sinfo, init_maps, retmap);
				num_dups = sum(node_dups.values());
				num_loss = sum(node_loss.values());
				gene_score = num_dups + num_loss;
				min_maps[gene_num].append([gene_score, num_dups, num_loss, maps, node_dups, node_loss]);
			else:
				gene_score = reconLCA(ginfo, sinfo, init_maps);

			total_score += gene_score;
			# Some counting.

		else:
			mt, minfo, hybrid_clade, hybrid_node, copy_node, = mul_tree[0], mul_tree[1], mul_tree[2], mul_tree[3], mul_tree[4];
			# Aggregate variables for the current GENE tree.

			gt_groups, gt_fixed = cur_groups[gene_num][0], cur_groups[gene_num][1];

			num_groups = len(gt_groups);
			# Retrieve gene tree info and collapsed groups for this gene tree-MUL-tree combo

			for combo in itertools.product(['','*'], repeat=num_groups):
			# We get all combinations of mappings for each node group. This is the time constraining step.
				group_map = [];
				for i in range(len(combo)):
					for node in gt_groups[i]:
						group_map.append(node + combo[i]);
				# This builds the current map for each group.

				for fixed in gt_fixed:
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

				if retmap:
					maps, node_dups, node_loss = reconLCA(ginfo, minfo, maps, retmap);
					num_dups = sum(node_dups.values());
					num_loss = sum(node_loss.values());
					cur_score = num_dups + num_loss;
					if cur_score <= gene_score:
						if cur_score < gene_score:
							gene_score = cur_score;
							min_maps[gene_num] = [];
						min_maps[gene_num].append([gene_score, num_dups, num_loss, maps, node_dups, node_loss])
				else:
					cur_score = reconLCA(ginfo, minfo, maps);
					if cur_score < gene_score:
						gene_score = cur_score;
				# Once the current maps have been initialized, we can simply call the normal LCA mapping algorithm
			## End mapping of one gene tree.

			total_score += gene_score;
	## End mapping all gene trees.

	if retmap:
		return min_maps;
	else:
		return mul_num, total_score;

# #############################################################################

# A couple ways to get the map combos:

# combo_ind = list(itertools.product(['','*'], repeat=len(node_ind)));
# if v == -2:
#	print "num combos", len(combo_ind);
# combos = list(itertools.product(['','*'], repeat=len(node_ind)));

# Old loading:
# if v == 0 and numiters > 100:
# 	numbars, donepercent = RC.loadingBar(itercount, numiters, donepercent, numbars);
# 	itercount = itercount + 1;
# # Only the loading bar displays when the program is running if -v is set to 0.






















