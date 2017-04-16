import sys, os, re

indir = sys.argv[1];
outfilename = sys.argv[2];

outfile = open(outfilename, "w");

filelist = os.listdir(indir);

for each in filelist:
	if each.find(".tre") == -1:
		continue;
	infilename = os.path.join(indir, each);
	infile = open(infilename, "r");
	intree = infile.read();
	infile.close();

	intree = re.sub('[)][\w]+;',');', intree);

	outfile.write(intree + "\n");

outfile.close();