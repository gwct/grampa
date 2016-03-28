#!/usr/bin/python
#############################################################################
# The main algorithmic functions for MUL-reconciliation mapping.
# Gregg Thomas
# Fall 2015, Combo algorithm implemented Spring 2016
#############################################################################

import sys, re, recontree as RT, itertools

#############################################################################

def reconLCA(lca_ginfo, sinfo, lca_maps):
# The LCA reconciliation mapping algorithm.

	dups = {};
	for g in lca_ginfo:
		dups[g] = 0;
	#The node classification of the gene tree. [key]:[value] -> [gene tree node]:[1,0]

	while [] in lca_maps.values():
	# To get the maps, a single post-order traversal of the tree is done.
		for g in lca_ginfo:
			if lca_ginfo[g][3] == 'tip' or lca_maps[g] != []:
				continue;
			# If the node is a tip, a map is already defined.

			desc = RT.getDesc(g, lca_ginfo)
			if lca_maps[desc[0]] == [] or lca_maps[desc[1]] == []:
				continue;
			# Otherwise, if both descendants haven't yet been mapped, continue.

			g_clade = RT.getClade(g, lca_ginfo);
			clade_maps = [];
			for g_tip in g_clade:
				clade_maps.append(lca_maps[g_tip][0]);
			# Get the species in the clade of the current node. Then get all
			# the possible maps from those species.

			lca_maps[g].append(RT.LCA(clade_maps,sinfo)[0]);

			if lca_maps[g][0] == lca_maps[desc[0]][0] or lca_maps[g][0] == lca_maps[desc[1]][0]:
				dups[g] = dups[g] + 1;
			#Now, if the map of g is identical to one of its descendants, it is a duplication node.

	dups, mdups, numdups = mulDupCount(lca_ginfo, sinfo, lca_maps);
	numloss = mulLossCount(lca_ginfo, sinfo, lca_maps, dups);
	# Call the functions to count the number of duplications and losses.

	return lca_maps;
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

def mulDupCount(dc_ginfo, dc_minfo, dc_maps):
# Given two trees (dictionaries) and a mapping between them, this function
# counts the number of duplication nodes.

	dup_count = 0;
	gt_dups = {};
	mt_dups = {};
	for m in dc_minfo:
		mt_dups[m] = 0;
	# The node classification of the species tree. [key]:[value] -> [species tree node]:[# of duplications along that branch]

	for g in dc_ginfo:
		if dc_ginfo[g][3] == 'tip':
			gt_dups[g] = 0;
			continue;
		desc = RT.getDesc(g, dc_ginfo);
		if dc_maps[g][0] == dc_maps[desc[0]][0] or dc_maps[g][0] == dc_maps[desc[1]][0]:
			gt_dups[g] = 1;
			mt_dups[dc_maps[g][0]] = mt_dups[dc_maps[g][0]] + 1;
			dup_count += 1;
		else:
			gt_dups[g] = 0;

	return gt_dups, mt_dups, dup_count;

#############################################################################

def mulLossCount(lc_ginfo, lc_minfo, lc_maps, lc_dups):
# Given two trees (dictionaries), a mapping between them, and the duplication nodes
# , this function counts the number of losses

	loss_count = 0;
	for g in lc_ginfo:
		if lc_ginfo[g][3] == 'root':
			continue;
		curanc = lc_ginfo[g][1];
		ancdepth = len(RT.nodeDepth(lc_maps[curanc][0],lc_minfo));
		gdepth = len(RT.nodeDepth(lc_maps[g][0],lc_minfo));

		glosses = 0;
		glosses = gdepth - ancdepth - 1;

		if lc_dups[curanc] != 0:
			glosses = glosses + 1;
		loss_count = loss_count + glosses;

	return loss_count

#############################################################################

def mulRecon(hybrid_clade, mt, minfo, gt, ginfo, v, check_nums):
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

	sisters = {};

	mul_hybrid_node = [n for n in minfo if set(RT.getClade(n, minfo)) == set(hybrid_clade)][0];
	copy_clade = [c + "*" for c in hybrid_clade];
	mul_copy_node = [n for n in minfo if set(RT.getClade(n, minfo)) == set(copy_clade)][0];
	# The copy clade is defined.

	hybrid_anc = minfo[mul_hybrid_node][1];
	copy_anc = minfo[mul_copy_node][1];

	sisters[''] = getSis(hybrid_anc, mul_hybrid_node, copy_clade, minfo);
	sisters['*'] = getSis(copy_anc, mul_copy_node, hybrid_clade, minfo);
	# These lines get any sister species from the hybrid and copy clades in the MUL-tree and that
	# clade's corresponding map. If there are no sisters, it stores an empty list.

	if v == -2:
		print "sisters", sisters;

	tmp_singles = [];
	tmp_groups = {};
	tmp_fixed = {};
	# The three possible classification of nodes.

	for g in ginfo:
		if ginfo[g][3] == 'tip':
			if g[g.index("_")+1:] in hybrid_clade:
				tmp_singles.append(g);
	# First, all hybrid species nodes in the gene tree are added to the singles list.
	## GET SINGLETONS

	for g in ginfo:
		if ginfo[g][3] != 'tip':
		# Next, for any non-tip node, we find out if the species that define it can be grouped
		# or fixed.
			#print g;
			d1, d2 = RT.getDesc(g, ginfo);
			d1_clade = RT.getClade(d1, ginfo);
			d1_spec_clade = [spec[spec.index("_")+1:] for spec in d1_clade];
			d2_clade = RT.getClade(d2,ginfo);
			d2_spec_clade = [spec[spec.index("_")+1:] for spec in d2_clade];
			# The clades for the descendants of both nodes are retrieved, and their corresponding
			# species are stored.


			if all(s in hybrid_clade for s in d1_spec_clade) and all(s in hybrid_clade for s in d2_spec_clade):
			# If the descendants from both nodes are all hybrid clade species, then we may be able to group them.
				if not any(s in d2_spec_clade for s in d1_spec_clade):
				# However, only if there is not more than one copy of a species among the clades can they be grouped.
					cur_clade = RT.getClade(g, ginfo);
					tmp_groups[g] = cur_clade;

					for s in cur_clade:
						if s in tmp_singles:
							tmp_singles.remove(s);
					# If these nodes are grouped, we must remove them from the singles list.
			## CHECK GROUPINGS


			if all(s in hybrid_clade for s in d1_spec_clade) and all(d1_spec_clade.count(s) <= 1 for s in hybrid_clade):
				for poss_map in sisters:
					if all(s in sisters[poss_map] for s in d2_spec_clade):
						tmp_fixed[d1] = [d1_clade, poss_map];


			elif all(s in hybrid_clade for s in d2_spec_clade) and all(d2_spec_clade.count(s) <= 1 for s in hybrid_clade):
				for poss_map in sisters:
					if all(s in sisters[poss_map] for s in d1_spec_clade):
						tmp_fixed[d2] = [d2_clade, poss_map];
			# For either descendant clade, if it contains only hybrid species (each with a count of 1 or 0) and the other
			# descendant contains only sister species that match those from the hybrid or copy clade in the MUL-tree, then
			# this clade can be fixed with the corresponding map.
			## FIX CLADES

	##############

	groups = tmp_groups.copy();
	for group in tmp_groups:
		cur_nodes = RT.getCladeNode(group, ginfo);
		for node in cur_nodes:
			if node in groups:
				del groups[node];
	# Removes any redundant descendant groups in the groups list.

	for node in tmp_fixed:
		if node in groups:
			del groups[node];
	# Removes any nodes from the groups list that are fixed.

	if v == -2:
		print "singles:", tmp_singles;
		print "tmp groups:", tmp_groups;
		print "groups:", groups;
		print "fixed:", tmp_fixed;

	init_num_groups = len(groups);
	hc_nodes = [];
	for node in groups:
		hc_nodes.append(groups[node]);
	del groups;
	for single in tmp_singles:
		hc_nodes.append([single]);
	# Converts the groups structure from a dictionary to a list. The ancestral nodes aren't important so we can discard
	# that information.

	if v == -2:
		print 'nodes:', hc_nodes;

	###############

	node_ind = range(len(hc_nodes));
	num_groups = len(hc_nodes);
	if v == -2:
		print "num nodes", len(node_ind);

	if check_nums:
		return len(tmp_singles), init_num_groups, len(tmp_fixed), num_groups;
	# If --checknums is True, we don't do any of the hard calcs, we just return some numbers as info.

	if num_groups > 15:
		return 0, 0, 0

	#combo_ind = list(itertools.product(['','*'], repeat=len(node_ind)));
	#if v == -2:
	#	print "num combos", len(combo_ind);
	# We get all combinations of mappings for each node group. This is the time constraining step.

	for combo in itertools.product(['','*'], repeat=len(node_ind)):
		hc_map = [];
		for i in range(len(combo)):
			for node in hc_nodes[i]:
				hc_map.append(node + combo[i]);
		# This builds the current map for each group.

		for fixed in tmp_fixed:
			for node in tmp_fixed[fixed][0]:
				hc_map.append(node + tmp_fixed[fixed][1]);
		# This adds the fixed maps onto the current combination of group mappings.

		if v == -2:
			print "curmap:", hc_map;

		# Now we do LCA mapping for the current combination of maps for the hybrid clade species.
		maps = {};
		for g in ginfo:
			if ginfo[g][3] == 'tip':
				speclabel = g[g.index("_")+1:];
				if g in hc_map:
					maps[g] = [speclabel];
				# If the node is in a hybrid clade, use the map in the current combination.
				elif g + "*" in hc_map:
					maps[g] = [speclabel + "*"];
				else:
					maps[g] = [speclabel];
				# Otherwise, the map is just the species label.

			else:
				maps[g] = [];
			# And if the node is not a tip, the map is empty.
		# Initialization of the mapping dictionary, using the current combination of hybrid clade maps. 
		# [key]:[value] -> [gene tree node]:[map]

		maps = reconLCA(ginfo, minfo, maps)
		# Once the current maps have been initialized, we can simply call the normal LCA mapping algorithm

		dups, mdups, numdups = mulDupCount(ginfo, minfo, maps);
		numloss = mulLossCount(ginfo, minfo, maps, dups);
		# Call the functions to count the number of duplications and losses.

		cur_score = numdups + numloss;

		if cur_score < min_score:
			min_score = cur_score;
			min_dups = numdups;
			min_loss = numloss;
			min_maps = maps;
		# If the current total score is smaller than the previous min, save it as the current min.

	return min_dups, min_loss, min_maps;
	# Return the maps with the minimal mutation score (dup + loss score).







