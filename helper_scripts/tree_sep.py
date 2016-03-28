import sys, os

infilename = sys.argv[1];
outdir = sys.argv[2];
listfilename = sys.argv[3];

listfile = open(listfilename, "w");


i = 1;
for line in open(infilename):
	outfilename = os.path.join(outdir, str(i)+".tre");
	outfile = open(outfilename, "w");
	outfile.write(line.replace("\n",""));
	outfile.close();
	listfile.write(outfilename + "\n");
	i += 1;


listfile.close();