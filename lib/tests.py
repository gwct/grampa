import sys, os, reconcore as RC

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

python_cmd = "python" + sys.argv[1];

grampath = os.path.dirname(__file__)[:-3];
grampath_script = os.path.join(grampath, "grampa.py");
grampath_s = os.path.join(grampath, "data", "manual_species_tree.tre");
grampath_g = os.path.join(grampath, "data", "manual_gene_trees.txt");
grampa_out = os.path.join(grampath, "data", "tests_out.txt");
outfilename = os.path.join(grampath, "data", "tests_log.txt");
outfile = open(outfilename, "w");
outfile.write("");
outfile.close();

tmpfile = start + "_tests.tmp";

RC.printWrite(outfilename, 1, "\nRUNNING GRAMPA TESTS");
RC.printWrite(outfilename, 1, start + "\n");

tests = ["labeltree", "multree", "checknum", "main"]
errors = {};
for t in tests:
	errors[t] = "";
# errors = { t : "" for t in tests }
numpass = 0;

RC.printWrite(outfilename, 1, "1: --labeltree test.........");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -v -1 --labeltree 2> " + tmpfile);
errors, numpass = catchErr("labeltree", errors, numpass, tmpfile);

RC.printWrite(outfilename, 1, "2: --multree test...........");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -g " + grampath_g + " -o " + grampa_out + " -v -1 --multree 2> " + tmpfile);
errors, numpass = catchErr("multree", errors, numpass, tmpfile);

RC.printWrite(outfilename, 1, "3: --checknum test..........");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -g " + grampath_g + " -o " + grampa_out + " -v -1 --checknum 2> " + tmpfile);
errors, numpass = catchErr("checknum", errors, numpass, tmpfile);

RC.printWrite(outfilename, 1, "4: MAIN test................");
os.system(python_cmd + " " + grampath_script + " -s " + grampath_s + " -g " + grampath_g + " -o " + grampa_out + " -v -1 2> " + tmpfile);
errors, numpass = catchErr("main", errors, numpass, tmpfile);

if numpass == 4:
	RC.printWrite(outfilename, 1, "\nDone! All tests pass!\n");
else:
	RC.printWrite(outfilename, 1, "\n" + str(4 - numpass) + " tests failed!");
	print("Check the tests_log.txt file in the test folder for more info!");
	RC.printWrite(outfilename, 1, "\n");
	outfile = open(outfilename, "a");
	for test in tests:
		if errors[test] != "":
			outfile.write(test + " failed with the following error:\n");
			outfile.write(errors[test] + "\n\n");
	outfile.close();

os.system("rm " + tmpfile);