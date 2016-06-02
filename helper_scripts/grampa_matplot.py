import sys

if len(sys.argv) != 3 or "-h" in sys.argv:
	print("\n# Usage: grampa_plot.py [input file] [output file]");
	print("# ---> [input file] must be a grampa output file.")
	print("# ---> [output file] will be an html file with your plot.\n")
	sys.exit();

infilename = sys.argv[1];
outfilename = sys.argv[2];

try:
	import numpy as np
	import matplotlib.pyplot as plt
except:
	sys.exit("Missing some of the required modules (numpy or matplotlib)")

score_dict = {};

for line in open(infilename):
	if line[0] == "#":
		continue;

	line = line.strip().split("\t");
	score_dict[line[1] + "-" + line[2]] = int(line[4]);

sorted_keys = sorted(score_dict, key=score_dict.get)
sorted_vals = [];

for key in sorted_keys:
	sorted_vals.append(score_dict[key]);

print(sorted_keys);
x = list(range(len(sorted_vals)));
print(sorted_vals);
plt.xticks(x, sorted_keys);
locs, labels = plt.xticks()
plt.setp(labels, rotation=90)
plt.scatter(x, sorted_vals);

plt.show();