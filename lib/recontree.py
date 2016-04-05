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

import re, sys
# re is used to replace internal labels in trees when necessary.
# I label internal nodes as '<#>'

#############################################################################

def buildMultree(h,p,tree,tree_info):
# This function builds a MUL-tree from a normally labeled species tree.
# Input:
#		h - a node that will have its subtree copied into the MUL-tree
#		p - a node at which the copied subtree will be placed
#		tree - the original species tree WITH internal nodes labeled
#		tree_info - the tree info dictionary (returned by treeparse) from the original tree

	if tree_info[h][3] == 'tip':
		hybrid = h;
	else:
		hybrid = getSubtree(h,tree);
		hybrid = re.sub('<[\d]+>','',hybrid);
	# Gets the subtree of the hybrid node

	if tree_info[p][3] == 'tip':
		copy = p;
	else:
		copy = getSubtree(p,tree);
		copy = re.sub('<[\d]+>','',copy);
	# Gets the subtree of the copy node.

	if (copy in hybrid and copy != hybrid) or tree_info[p][3] == 'root':
		return "NULL";
	# Copy nodes shouldn't be within the hybrid subtree and shouldn't be at the root.

	if tree_info[h][3] != 'tip':
		for node in tree_info:
			if node in hybrid and tree_info[node][3] == 'tip':
				tree = tree.replace(node, node+"_1");
				copy = copy.replace(node, node+"_1")
	# Some re-labeling if necessary.
	tree = re.sub('<[\d]+>','',tree);

	if tree_info[h][3] == 'tip':
			mul_clade = "(" + copy + "," + hybrid + "*)";
	else:
			mul_clade = "(" + copy + "," + hybrid + ")";

	mul_tree = tree.replace(copy,mul_clade);
	# Combines the clades and replaces the copy clade in the original tree to create the MUL-tree.

	mul_tree = re.sub('<[\d]+>','',mul_tree);
	if tree_info[h][3] != 'tip':
		for node in tree_info:
			if node in hybrid and tree_info[node][3] == 'tip':
				mul_tree = mul_tree.replace(node, node+"*");
		mul_tree = mul_tree.replace("*_1", "");
	# Some relabling of the hybrid species.

	return mul_tree;

#############################################################################

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

				for a in xrange(len(indcheck)):
					if indcheck[a] == -1:
						indcheck[a] = 10000;

				curbranch = bltree[d+1:min(indcheck)];
				return curbranch;
		d = d + 1;
	startind = d;

#############################################################################

def remBranchLength(treestring):
# Removes branch lengths from a tree.
	if "):" in treestring:
		treestring = re.sub(':[\d.Ee-]+[)]', ')', treestring);
	if ":" in treestring:
		treestring = re.sub('[)][\d.-]+:[\d.Ee-]+', ')', treestring);
		treestring = re.sub(':[\d.Ee-]+', '', treestring);
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
		if c_treedict[d][3] != 'tip':
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
		if c_treedict[d][3] != 'tip':
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
	while tree_dict[node][3] != 'root':
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
		if treedict[spec][3] == 'root':
			continue;
		curanc = treedict[spec][1];
		ancs[spec].append(curanc);
		while treedict[curanc][3] != 'root':
			curanc = treedict[curanc][1];
			ancs[spec].append(curanc);

	intersect_anc = set.intersection(*map(set, ancs.values()))
	lcp = [t for t in ancs.values()[0] if t in intersect_anc]

	#lcp = sorted(set.intersection(*map(set, ancs.values())), key=lambda x: ancs.values()[0].index(x))
	monophyletic = 0;
	if set(getClade(lcp[0],treedict)) == set(spec_list):
		monophyletic = 1;

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

def nodeDepth(n_spec, n_treedict):
# This function returns a list of the nodes between the current node and the
# root of the tree. Therefore, the node depth is the length of this list.
# Larger depths are more towards the tips of the tree. The root node has
# depth 0.

	if n_treedict[n_spec][3] == 'root':
		return [];

	ancs = []
	curanc = n_treedict[n_spec][1];
	ancs.append(curanc);
	while n_treedict[curanc][3] != 'root':
		curanc = n_treedict[curanc][1];
		ancs.append(curanc);
	return ancs;

#############################################################################

def treeParseNew(tree, tree_type):
# The treeParse function takes as input a rooted phylogenetic tree with branch lengths and returns the tree with node labels and a
# dictionary with usable info about the tree in the following format:
# New (current) format:
# node:[branch length, ancestral node, ancestral branch length, node type]
#
# Tree type 1: tree has branch lengths.
# Tree type 2: tree is just topology.

	tree = tree.replace("\n","");
	if tree[len(tree)-1:] != ";":
		tree = tree + ";";
	tree = re.sub('\)[\d\w.]*:','):',tree);
	## Some string handling

	new_tree = "";
	z = 0;
	numnodes = 1;

	while z < (len(tree)-1):
		if tree_type == 1:
			if tree[z] == ":" and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		if tree_type == 2:
			if (tree[z] == "," or tree[z] == ")") and tree[z-1] == ")":
				new_tree = new_tree + "<" + str(numnodes) + ">";
				numnodes = numnodes + 1;
		new_tree = new_tree + tree[z];
		z = z + 1;
	if new_tree[-1] == ")":
		rootnode = "<" + str(numnodes) + ">"
		new_tree = new_tree + rootnode;
	else:
		rootnode = new_tree[new_tree.rfind(")")+1:];
	## This first block labels all internal nodes with the format <#>
	#print "-----------------------------------";
	#print new_tree;
	#print "-----------------------------------";

	ancs = {};
	nofo = {};

	z = 0;
	startind = 0;
	while z < (len(new_tree)-1):
	##Here, the ancestral nodes of each node are found

		if tree_type == 1:
		## The major difference between trees with branch lengths (type 1) and without (type 2) is seen here. Finding the ancestral nodes requires
		## almost a completely different set of logic statements.
			if new_tree[z] == ":":
				curnode = new_tree[max(new_tree.rfind("(",startind,z),new_tree.rfind(")",startind,z),new_tree.rfind(",",startind,z))+1:z];
				numcpneeded = 1
				numcp = 0;
				nofo[curnode] = [];

				a = z;

				while a < (len(new_tree)-1):
					if new_tree[a] == "(":
						numcpneeded = numcpneeded + 1;
					if new_tree[a] == ")" and numcpneeded != numcp:
						numcp = numcp + 1;
					if new_tree[a] == ")" and numcpneeded == numcp:
						#if a == (len(new_tree)-5):
						#	curanc = new_tree[a+1:];
						if new_tree[a+1:].find(":") == -1:
							#curanc = new_tree[len(new_tree)-5:];
							curanc = new_tree[new_tree.rfind(")")+1:]
						else:
							curanc = new_tree[a+1:new_tree.index(":", a)];
						a = 100000000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		if tree_type == 2:
			if new_tree[z] == "," or new_tree[z] == ")":
				curnode = new_tree[max(new_tree.rfind("(",startind,z),new_tree.rfind(")",startind,z),new_tree.rfind(",",startind,z))+1:z];
				numcpneeded = 1
				numcp = 0;
				nofo[curnode] = [];

				a = z;

				while a < (len(new_tree)-1):
					if new_tree[a] == "(":
						numcpneeded = numcpneeded + 1;
					if new_tree[a] == ")" and numcpneeded != numcp:
						numcp = numcp + 1;
					if new_tree[a] == ")" and numcpneeded == numcp:
						if a == (len(new_tree)-4):
							curanc = new_tree[a+1:];
						else:
							mindex = 999999999;
							for c in ["(",")",","]:
								cind = new_tree.find(c,a+1);
								if cind < mindex and cind != -1:
									mindex = cind;
									minchar = c;
							curanc = new_tree[a+1:mindex];
						a = 10000;

						ancs[curnode] = curanc;
					a = a + 1;
				startind = z;

		z = z + 1;
	## End ancestral node block

	#for key in ancs:
	#	print key + ":", ancs[key]
	#print "---------";
	#sys.exit()

	## The next block gets all the other info for each node. This is easy now that the ancestral nodes are stored.
	nofo[rootnode] = [];
	for node in nofo:
		if tree_type == 1:
			cur_bl = getBranchLength(new_tree,node);
		elif tree_type == 2:
			if node == rootnode:
				cur_bl = None;
			else:
				cur_bl = "NA";
		nofo[node].append(cur_bl);

		if node != rootnode:
			cur_anc = ancs[node];
			nofo[node].append(cur_anc);
			if tree_type == 1:
				cur_anc_bl = getBranchLength(new_tree,cur_anc);
			elif tree_type == 2:
				cur_anc_bl = "NA";
			nofo[node].append(cur_anc_bl);
		else:
			j = 0;
			while j < 2:
				nofo[node].append("");
				j = j + 1;

	for node in nofo:
		if node == rootnode:
			nofo[node].append("root");
		elif getDesc(node,nofo) == [node]:
			nofo[node].append("tip");
		else:
			nofo[node].append("internal");

	## End info retrieval block.

#	for key in nofo:
#		print key + ":" + str(nofo[key]);


	return nofo, new_tree;

#############################################################################