import sys
import os
import subprocess
import grampa_lib.reconcore as RC

###############################

def catchErr(test, err_dict, p, t_file):
	err = open(t_file, "r").read();
	if err != "":
		errors[test] = err;
		RC.printWrite(outfilename, 1, "Failed!");
	else:
		RC.printWrite(outfilename, 1, "OK");
		p += 1;
	return err_dict, p;

###############################

start = RC.getLogTime();

python_cmd = sys.argv[1];
grampath = os.path.dirname(__file__)[:-10];
grampath_script = os.path.join(grampath, "grampa.py");
grampath_s = os.path.join(grampath, "data", "manual_species_tree.tre");
grampath_g = os.path.join(grampath, "data", "manual_gene_trees.txt");
grampa_out = os.path.join(grampath, "tests_out_3akjg4z");
outfilename = os.path.join(grampath, "tests_log_" + start + ".txt");
#RC.filePrep(outfilename);

with open(outfilename, "w") as outfile:
	outfile.write("");

tmpfile = "tests_" + start + ".tmp";

RC.printWrite(outfilename, 2, "\nRUNNING GRAMPA TESTS");
RC.printWrite(outfilename, 2, start + "\n");

tests = ["labeltree", "multree", "checknum", "main"]
errors = {};
for t in tests:
	errors[t] = "";
numpass = 0;

RC.printWrite(outfilename, 2, "1: --labeltree test.........");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -v 0 --labeltree 2> " + tmpfile);
errors, numpass = catchErr("labeltree", errors, numpass, tmpfile);
os.system("rm -rf " + grampa_out);

RC.printWrite(outfilename, 2, "2: --buildmultrees test......");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -g " + grampath_g + " -o " + grampa_out + " -v 0 --buildmultrees 2> " + tmpfile);
errors, numpass = catchErr("multree", errors, numpass, tmpfile);
os.system("rm -rf " + grampa_out);

RC.printWrite(outfilename, 2, "3: --checknum test..........");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -g " + grampath_g + " -o " + grampa_out + " -v 0 --checknums 2> " + tmpfile);
errors, numpass = catchErr("checknum", errors, numpass, tmpfile);
os.system("rm -rf " + grampa_out);

RC.printWrite(outfilename, 2, "4: MAIN test................");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -g " + grampath_g + " -o " + grampa_out + " -v 0 --maps 2> " + tmpfile);
errors, numpass = catchErr("main", errors, numpass, tmpfile);
os.system("rm -rf " + grampa_out);

if numpass == 4:
	RC.printWrite(outfilename, 2, "\nDone! All tests pass!\n");
else:
	RC.printWrite(outfilename, 2, "\n" + str(4 - numpass) + " tests failed!");
	print("Check the " + outfilename + " file for more info!");
	RC.printWrite(outfilename, 2, "\n");
	outfile = open(outfilename, "a");
	for test in tests:
		if errors[test] != "":
			outfile.write(test + " failed with the following error:\n");
			outfile.write(errors[test] + "\n\n");
	outfile.close();

os.system("rm " + tmpfile);