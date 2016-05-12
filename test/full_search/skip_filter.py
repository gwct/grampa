import sys

infilename = "h1_20_det.txt"

trees = range(1,4624);
to_remove = [];
i = 1;

print "# Getting trees to remove...";
for line in open(infilename):
	if line[:2] != "GT":
		continue;

	if "Skipping" in line:
		gt = line[line.index("-")+1:line.index(" ")];
		# print line;
		# print gt;
		# print "---"
		if int(gt) in trees:
			trees.remove(int(gt));
		if int(gt) not in to_remove:
			to_remove.append(int(gt));
print "#", len(to_remove), "trees to remove.";

treefilename = "../yeast206_trees_rfbr.txt";
outfilename = "yeast206_skipfilter.txt";

print "# Reading trees...";
treelist = open(treefilename, "r").readlines();
print "#", len(treelist), "trees read.";

print "# Filtering trees...";
newlist = [];
for t in range(len(treelist)):
	if t not in to_remove:
		newlist.append(treelist[t]);

print "# Writing remaining", len(newlist), "trees to output file...";
outfile = open(outfilename, "w");
for tree in newlist:
	outfile.write(tree);
outfile.close();

print "# Done!";