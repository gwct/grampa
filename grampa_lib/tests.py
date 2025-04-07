import sys
import os
import subprocess
import shutil
import grampa_lib.reconcore as RC

def run_test(args, output_dir):
    """
    Run a test using subprocess and capture stderr output.
    Returns an empty string if no error occurred, or the error message if one did.
    Also removes the output directory if it exists.
    """
    try:
        print(f"Running command: {' '.join(args)}")
        result = subprocess.run(args, capture_output=True, text=True, check=False);
        error_message = result.stderr.strip();
        #print(result.stderr)
    except Exception as e:
        error_message = str(e);
    # Run the command and capture the output
    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir);
    # Remove the output directory after the test command
    
    return error_message;

def main(test_opt):
    if test_opt == "bioconda":
        python_cmd = "";
        grampath = ".";
        grampath_script = "grampa";
        test_dir = "bioconda-test-data";
    # biconda tests
    else:
        python_cmd = sys.argv[1];
        grampath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."));
        grampath_script = os.path.join(grampath, "grampa.py");
        # Get the absolute path of the grampa.py script, which is located in the parent directory of the current script

        test_dir = os.path.join(grampath, "data", "bioconda-test-data");
        # Get the path to the test directory, which is located in the data subdirectory of the grampa directory
    # Local tests

    if not os.path.exists(test_dir):
        print(f"Error: Test directory '{test_dir}' not found.", file=sys.stderr);
        sys.exit(1);

    grampath_s = os.path.join(test_dir, "manual_species_tree.tre");
    grampath_g = os.path.join(test_dir, "manual_gene_trees.txt");
    # Get the paths to the species tree and gene trees files

    grampa_out = os.path.join(grampath, "grampa-tests-out-3akjg4z");
    # Create a temporary output directory for the tests
    
    print("\nRUNNING GRAMPA TESTS\n");
    
    tests = [
        {
            "display": "--labeltree test",
            "args": [python_cmd, grampath_script, "-s", grampath_s, "-v", "0", "--labeltree"]
        },
        {
            "display": "--buildmultrees test",
            "args": [python_cmd, grampath_script, "-s", grampath_s, "-g", grampath_g,
                    "-o", grampa_out, "-v", "0", "--buildmultrees"]
        },
        {
            "display": "--checknum test",
            "args": [python_cmd, grampath_script, "-s", grampath_s, "-g", grampath_g,
                    "-o", grampa_out, "-v", "0", "--checknums"]
        },
        {
            "display": "main test",
            "args": [python_cmd, grampath_script, "-s", grampath_s, "-g", grampath_g,
                    "-o", grampa_out, "-v", "0", "--maps"]
        }
    ];
    # List of tests to run, each with a display name and command line arguments
    # The display name is used for printing the test results
    # The command line arguments are passed to the subprocess call
    # The output directory is created and removed for each test
    
    if test_opt == "bioconda":
        for i in range(len(tests)):
            tests[i]["args"].pop(0);
    # Remove the first argument (python_cmd) from the command line arguments
    # This is necessary for the bioconda tests, where the script is run directly

    errors = {};
    numpass = 0;
    fixed_width = 30;  # Adjust this width to get the desired number of dots
    
    for idx, test in enumerate(tests, start=1):
        test_desc = f"{idx}: {test['display']}"
        padded_desc = f"{test_desc:.<{fixed_width}}"
        # Print the description immediately without newline
        sys.stdout.write(padded_desc)
        sys.stdout.flush()
        # Run the test
        error_message = run_test(test["args"], grampa_out)
        result_str = "OK" if not error_message else "FAILED"
        # Print the result on the same line
        print(result_str)
        errors[test['display']] = error_message
        if not error_message:
            numpass += 1

    if numpass != len(tests):
        print("\nSome tests failed. Please review the error details below:\n")
        for test in tests:
            if errors[test["display"]]:
                print(f"{test['display']} failed with the following error:")
                print(errors[test["display"]])
                print()
                print("+++")
        sys.exit(1);
    # If any test failed, print the error messages and exit with a non-zero status
    else:
        print("\nDone! All tests pass!\n")
    # Otherwise, print a success message

if __name__ == '__main__':
    test_opt = "";
    if "--bioconda" in sys.argv:
        test_opt = "bioconda";
    main(test_opt);
