import sys

if len(sys.argv) != 3 or "-h" in sys.argv:
	print("\n# This is a beta version of this script and may be buggy.")
	print("# Usage: grampa_plot.py [input file] [output file]");
	print("# ---> [input file] must be a grampa output (_out.txt) file.")
	print("# ---> [output file] will be a png image with your plot.\n")
	sys.exit();

infilename = sys.argv[1];
outfilename = sys.argv[2];
if outfilename[len(outfilename)-4:] != ".png":
	outfilename += ".png";

try:
	import numpy as np
	import matplotlib.pyplot as plt
except:
	sys.exit("Missing some of the required modules (numpy or matplotlib)")

score_dict = {};

for line in open(infilename):
	if line[0] == "#" or "The" in line or "Score" in line:
		continue;

	line = line.strip().split("\t");
	if len(line) == 4:
		score_dict[line[0]] = int(line[3]);
	else:
		score_dict[line[1] + "-" + line[2]] = int(line[4]);

sorted_keys = sorted(score_dict, key=score_dict.get)
sorted_vals = [];

for key in sorted_keys:
	sorted_vals.append(score_dict[key]);

x = list(range(len(sorted_vals)));
fig = plt.figure();
fig.set_size_inches(12, 6)
plt.xticks(x, sorted_keys);
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)
plt.scatter(x, sorted_vals);
fig.suptitle('GRAMPA Results: ' + infilename);
plt.xlabel('H1-H2 Node');
plt.ylabel('Score');
fig.savefig(outfilename, bbox_inches='tight');
plt.show();