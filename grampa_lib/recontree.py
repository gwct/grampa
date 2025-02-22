#!/usr/bin/python
#############################################################################
# Functions for retrieving information from newick formatted, rooted phylogenetic trees.
# Additional functions added for use with reconciliation with MUL-trees.
#
# The core function here is treeParse. It takes a newick-formatted tree as a string
# as input and returns the same tree with internal nodes labeled and a dictionary of information
# about each node. It also takes a tree type input parameter, with 1 being a tree with branch
# lengths and 2 being a tree without branch lengths.
#
# The [key] : [value] format of the tree dictionary is:
#	[node] : [branch length, ancestral node, ancestral branch length, node type]
#
# This tree dictionary is used by almost all other tree functions.
#
# Gregg Thomas
# Spring 2013-present
# Forked from core treeparse 12.08.2015
#############################################################################

import re
import sys
# re is used to replace internal labels in trees when necessary.
# I label internal nodes as '<#>'

def getBranchLength(bltree, spec_label):
# Returns the branch length of a species given a newick formatted tree. Used by treeParse.
	d = 0;
	startind = 0;

	while d < (len(bltree)-1):
		if bltree[d] == ":":
			current_node = bltree[max(bltree.rfind("(",startind,d),bltree.rfind(")",startind,d),bltree.rfind(",",startind,d))+1:d];
			if current_node == spec_label:

				opind = bltree.find("(",d);
				cpind = bltree.find(")",d);
				coind = bltree.find(",",d);

				indcheck = [opind,cpind,coind];

				for a in range(len(indcheck)):
					if indcheck[a] == -1:
						indcheck[a] = 10000;

				curbranch = bltree[d+1:min(indcheck)];
				return curbranch;
		d = d + 1;
	startind = d;

#############################################################################

def remBranchLength(treestring):
# Removes branch lengths from a tree.

	treestring = re.sub('[)][\d\w<>.eE_:-]+', ')', treestring);
	treestring = re.sub(':[\d.eE-]+', '', treestring);
	# treestring = re.sub('[)][\d.eE-]+:[\d.eE-]+', ')', treestring);
	# treestring = re.sub(':[\d.eE-]+', '', treestring);
	# treestring = re.sub('<[\d]+>', '', treestring);

	return treestring;

#############################################################################

def getDesc(d_spec, d_treedict):
# This function takes a species in the current tree and the dictionary of the current tree
# (returned by treeparse) and finds the direct descendants of the species.

	d_list = [];
	for node in d_treedict:
		if d_treedict[node][1] == d_spec:
			d_list.append(node);

	if d_list == []:
		return [d_spec];
	# If no descendants have been found, the node must be a tip and the function
	# returns the same node.
	else:
		return d_list;

#############################################################################

def getClade(c_spec, c_treedict):
# This function takes a node in the current tree and the dictionary of the current tree
# (returned by treeparse) and finds all tip labels that are descendants of the current node.
# This is done by getting the direct descendants of the current node with getDesc and then
# recursively calling itself on those descendants.

	clade = [];
	c_desc = getDesc(c_spec, c_treedict);
	for d in c_desc:
		if c_treedict[d][2] != 'tip':
			clade.append(getClade(d, c_treedict));
			# Recursion
		else:
			clade.append(d);

	r_clade = [];
	for c in clade:
		if type(c) == list:
			for cc in c:
				r_clade.append(cc);
		else:
			r_clade.append(c);

	return r_clade;

#############################################################################

def getCladeNode(c_spec, c_treedict):
# This function takes a node in the current tree and the dictionary of the current tree
# (returned by treeparse) and finds all tip labels that are descendants of the current node.
# This is done by getting the direct descendants of the current node with getDesc and then
# recursively calling itself on those descendants.

	clade = [];
	c_desc = getDesc(c_spec, c_treedict);
	for d in c_desc:
		if c_treedict[d][2] != 'tip':
			clade.append(getCladeNode(d, c_treedict));
			# Recursion
		clade.append(d);

	r_clade = [];
	for c in clade:
		if type(c) == list:
			for cc in c:
				r_clade.append(cc);
		else:
			r_clade.append(c);

	return r_clade;

#############################################################################

def pathToRoot(node, tree_dict):
	ptr = [node];
	while tree_dict[node][2] != 'root':
		ptr.append(tree_dict[node][1]);
		node = tree_dict[node][1];
	return ptr;

#############################################################################

def LCA(spec_list, treedict):
# LCA stands for Least Common Ancestor, and this function takes a list of species
# and a tree dictionary (returned by treeParse) and finds the first node that has 
# all of those species as descendants (the LCA).
# This is accomplished by storing a list of all nodes from each node to the root.
# Then, the common path from root to a node is found and the first element of this
# is returned as the LCA.
# The function also checks if the input species list is monophyletic by checking
# if the clade of the LCA is the same as the input list.

	ancs = {};
	for spec in spec_list:
		ancs[spec] = [spec];

	for spec in spec_list:
		if treedict[spec][2] == 'root':
			continue;
		curanc = treedict[spec][1];
		ancs[spec].append(curanc);
		while treedict[curanc][2] != 'root':
			curanc = treedict[curanc][1];
			ancs[spec].append(curanc);

	intersect_anc = set.intersection(*list(map(set, list(ancs.values()))))
	lcp = [t for t in list(ancs.values())[0] if t in intersect_anc]

	#lcp = sorted(set.intersection(*map(set, ancs.values())), key=lambda x: ancs.values()[0].index(x))
	monophyletic = False;
	if set(getClade(lcp[0],treedict)) == set(spec_list):
		monophyletic = True;

	return lcp[0], monophyletic;

#############################################################################

def getSubtree(node, tree):
# Given a node and a tree with internal nodes labeled, this returns the
# subtree rooted at that node.

	subtree = "";
	partree = tree[:tree.index(node)][::-1];
	cp = 0;
	op = 0;
	for c in partree:
		if c == ")":
			cp = cp + 1;
		if c == "(":
			op = op + 1;
		subtree = subtree + c;
		if cp == op:
			break;
	return subtree[::-1];

#############################################################################

def getSis(s_node, s_dict):
# Gets the sister node of a given node. Sister node is defined as the other
# descendent node of their direct common ancestor node.

	s_anc = s_dict[s_node][1];
	anc_desc = getDesc(s_anc, s_dict);
	sis_node = [node for node in anc_desc if node != s_node];
	return "".join(sis_node);

#############################################################################

def nodeDepth(n_spec, n_treedict):
# This function returns a list of the nodes between the current node and the
# root of the tree. Therefore, the node depth is the length of this list.
# Larger depths are more towards the tips of the tree. The root node has
# depth 0.

	if n_treedict[n_spec][2] == 'root':
		return [];

	ancs = []
	curanc = n_treedict[n_spec][1];
	ancs.append(curanc);
	while n_treedict[curanc][2] != 'root':
		curanc = n_treedict[curanc][1];
		ancs.append(curanc);
	return ancs;

#############################################################################

def sortNodes(tree_dict):
	sorted_nodes = [];
	for n in tree_dict:
		if tree_dict[n][2] != 'tip':
			sorted_nodes.append(int(n.replace("<","").replace(">","")));
	sorted_nodes = sorted(sorted_nodes);

	return sorted_nodes;

#############################################################################

def treeParse(tree, debug=0):
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# New (current) format:
# node:[branch length (if present), ancestral node, node type, node label (if present)]

	tree = tree.strip();
	if tree[-1] != ";":
		tree += ";";
	# Some string handling

	nodes = {};
	bl = {};
	supports = {};
	ancs = {};
	# Initialization of all the tracker dicts

	topology = remBranchLength(tree);
	nodes = {};
	for n in topology.replace("(","").replace(")","").replace(";","").split(","):
		nodes[n] = 'tip';
	# nodes = { n : 'tip' for n in topology.replace("(","").replace(")","").replace(";","").split(",") };
	# Retrieval of the tip labels

	new_tree = "";
	z = 0;
	numnodes = 1;
	while z < (len(tree)-1):
		new_tree += tree[z];
		if tree[z] == ")":
			node_label = "<" + str(numnodes) + ">";
			new_tree += node_label;
			nodes[node_label] = 'internal';
			numnodes += 1;
		z += 1;
	nodes[node_label] = 'root';
	# This labels the original tree as new_tree and stores the nodes and their types in the nodes dict

	topo = "";
	z = 0;
	numnodes = 1;
	while z < (len(topology)-1):
		topo += topology[z];
		if topology[z] == ")":
			node_label = "<" + str(numnodes) + ">";
			topo += node_label;
			numnodes += 1;
		z += 1;
	# This labels the topology with the same internal labels

	for node in nodes:
	# One loop through the nodes to retrieve all other info

		if node + ")" in tree or node + "," in new_tree:
		# If the node is followed immediately by a ) or , then there are no branch lengths of supports to collect
			supports[node] = "NA";
			bl[node] = "NA";

		elif node + ":" in new_tree:
		# If the node is followed immediately by a : then there is a branch length, but no support, to collect
			supports[node] = "NA";
			cur_bl = re.findall(node + ":[\d.Ee-]+", new_tree);
			cur_bl = cur_bl[0].replace(node + ":", "");
			bl[node] = cur_bl;

		else:
		# Otherwise we must collect both support and branch length or just support
			cur_bsl = re.findall(node + "[\d*+.Ee-]+:[\d.Ee-]+", new_tree);
			if cur_bsl:
			# If the pattern above is found then the node has both support and branch length
				cur_bs = cur_bsl[0].replace(node, "");
				cur_bs = cur_bs[:cur_bs.index(":")];
				cur_bl = cur_bsl[0].replace(node + cur_bs + ":", "");
				supports[node] = cur_bs;
				bl[node] = cur_bl;
			else:
			# If it is not found then the branch only has a support label
				cur_bs = re.findall(node + "[\w*+.<> -]+", new_tree);
				supports[node] = cur_bs;
				bl[node] = "NA";

		# Next we get the ancestral nodes. If the node is the root this is set to NA.
		if nodes[node] == 'root':
			ancs[node] = "NA";
			continue;

		match_node = node;
		if "*" in node:
			match_node = match_node.replace("*", "\*");
		if "^" in node:
			match_node = match_node.replace("^", "\^");
		if "+" in node:
			match_node = match_node.replace("+", "\+");

		anc_match = re.findall('[(),]' + match_node + '[(),]', new_tree);
		anc_tree = new_tree[new_tree.index(anc_match[0]):][1:];
		if anc_match[0][-1] != "*" and anc_tree[len(anc_match[0])] == "*":
			anc_tree = new_tree[new_tree.rindex(anc_match[0]):][1:];
		# Ancestral labels are always to the right of the node label in the text of the tree, so we start our scan from the node label

		if debug == 1:
			print(node);
			print(anc_match);
			print(anc_tree);
			print("---");

		cpar_count = 0;
		cpar_need = 1;

		for i in range(len(anc_tree)):
		# We find the ancestral label by finding the ) which matches the nesting of the number of ('s found
			if anc_tree[i] == "(":
				cpar_need = cpar_need + 1;
			if anc_tree[i] == ")" and cpar_need != cpar_count:
				cpar_count = cpar_count + 1;
			if anc_tree[i] == ")" and cpar_need == cpar_count:
				anc_tree = anc_tree[i+1:];
				ancs[node] = anc_tree[:anc_tree.index(">")+1];
				break;

	nofo = {};
	for node in nodes:
		if nodes[node] == "tip":
			orig_node = node;
			while node[0] == " " or node[-1] == " ":
				if node[0] == " ":
					node = node[1:];
				if node[-1] == " ":
					node = node[:-1];
			nofo[node] = [bl[orig_node], ancs[orig_node], nodes[orig_node], supports[orig_node]];
			topo = topo.replace(orig_node, node);
		# This block checks for preceding and trailing spaces in the tip labels of the tree and removes them
		# if found. This was causing errors later on when matching gene tree to species tree labels for mapping.

		else:
			nofo[node] = [bl[node], ancs[node], nodes[node], supports[node]];
		# For all internal nodes just add the node as is.
	# Now we just restructure everything to the old format for legacy support

	if debug == 1:
	# Debugging options to print things out
		print("\ntree:\n" + tree + "\n");
		print("new_tree:\n" + new_tree + "\n");
		print("topology:\n" + topo + "\n");
		print("nodes:");
		print(nodes);
		print()
		print("bl:");
		print(bl);
		print()
		print("supports:");
		print(supports);
		print()
		print("ancs:");
		print(ancs);
		print()
		print("-----------------------------------");
		print()
		print("nofo:");
		print(nofo);
		print("")

	return nofo, topo;

#############################################################################
