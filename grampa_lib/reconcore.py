#############################################################################
# Reconciliation CORE functions
# Gregg Thomas
# August 2013-present
# Forked from core on 12.13.2015
#############################################################################

import sys
import os
import subprocess
import datetime
import time
import timeit

#############################################################################

def spacedOut(string, totlen, sep=" "):
# Properly adds spaces to the end of a message to make it a given length
    spaces = sep * (totlen - len(string));
    return string + spaces;

#############################################################################

def isPosInt(numstr, default=False, minval=1, maxval=False):
# Check if a string is a positive integer
    try:
        num = int(numstr);
    except:
        return default;

    if num < minval:
        return default;
    elif maxval and num > maxval:
        return default;
    else:
        return num;

#############################################################################

def getOutTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m-%d-%Y.%I-%M-%S");

#############################################################################

def getLogTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%I.%M.%S");
    # return datetime.datetime.now().strftime("%m.%d.%Y-%I.%M.%S");

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y  %H:%M:%S");

#############################################################################

def getDate():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%m.%d.%Y");

#############################################################################

def getTime():
# Function to get the date and time in a certain format.
    return datetime.datetime.now().strftime("%H:%M:%S");

#############################################################################

def printWrite(o_name, v, o_line1, o_line2="", pad=0):
#Function to print a string AND write it to the file.
    if o_line2 == "":
        outline = o_line1;
    else:
        outline = o_line1 + "."*(pad-len(o_line1)) + o_line2;
    if v in [2,3]:
        print(outline);
    if v != -1:
        f = open(o_name, "a");
        f.write(outline + "\n");
        f.close();

#############################################################################

def report_step(globs, step, step_start_time, step_status, start=False, full_update=False, out_v=1):
# Uses psutil to gather memory and time info between steps and print them to the screen.

    dashes = 150
    if globs['psutil']:
        import psutil;
        dashes = 175;
    # Determine the number of dashes to frame the update table depending on the presence of psutil

    cur_time = timeit.default_timer();
    # The time at the start of the status update

    if globs['num-opt']:
        return cur_time;

    col_widths = [ 14, 10, 50, 40, 20, 16 ];
    if globs['psutil']:
        col_widths += [18, 10];
    # The column widths

    if start:
        headers = [ "# Date", "Time", "Current step", "Status", "Elapsed time (s)", "Step time (s)" ];
        if globs['psutil']:
            headers += ["Current mem (MB)", "Virtual mem (MB)"]
        # A list of the headers

        headers = "".join([ spacedOut(str(headers[i]), col_widths[i]) for i in range(len(headers)) ]);
        # Converting the list to a string based on the column widths

        printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
        printWrite(globs['logfilename'], globs['log-v'], headers);
        printWrite(globs['logfilename'], globs['log-v'], "# " + "-" * dashes);
        # Print the dashes and the headers
    # The first call is just to print the headers

    ##########

    else:
        prog_elapsed = str(round(cur_time - globs['starttime'], 5));
        # Get the total amount of time that the program has been running

        if not step_start_time:
        # If no step start time is given, then this is the first entry for this status
        # update, that will display "In progress..." or similar.

            out_line = [ "# " + getDate(), getTime(), step, step_status ];
            # The output for the initial status entry includes the date, time, step label, and progress message

            term_col_widths = col_widths[:4];
            # Only get the first 4 column widths for the initial status entry

            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];

            if full_update:
                out_line += "\n";
            # For some status updates, intermediate info will be printed, in which case we add a newline here

            if not globs['log-v'] <= 1: 
                sys.stdout.write("".join(out_line));
                sys.stdout.flush();
            # Convert the output list to a string, write, and flush stdout

        ## The initial status entry to display "In progress..."
        #####

        else:
            step_elapsed = str(round(cur_time - step_start_time, 5));
            # Get the full step time here

            out_line = [ step_status, prog_elapsed, step_elapsed ];
            # Gather info for the full output line to print to screen

            if globs['psutil']:
                mem = round(sum([p.memory_info()[0] for p in globs['pids']]) / float(2 ** 20), 5);
                vmem = round(sum([p.memory_info()[1] for p in globs['pids']]) / float(2 ** 20), 5);
                out_line += [str(mem), str(vmem)];
            # If psutil is present, get current memory info

            term_col_widths = col_widths[3:];
            # Get the column widths for the print to screen output

            file_line = [ "# " + getDate(), getTime(), step ] + out_line;
            file_col_widths = col_widths[:3] + [30] + col_widths[4:];
            # For output to the file, we write the whole line each time
            # Add the initial entry fields here
            # This will also be used for some status updates where the whole message needs to be printed
            # to the screen
            
            out_line = [ spacedOut(str(out_line[i]), term_col_widths[i]) for i in range(len(out_line)) ];
            file_line = [ spacedOut(str(file_line[i]), col_widths[i]) for i in range(len(file_line)) ];
            # Compile both the truncated and the full status update

            if not globs['log-v'] <= 1:
                if full_update:
                    sys.stdout.write("".join(file_line) + "\n");
                    sys.stdout.flush();
                else:         
                    sys.stdout.write("\b" * 40);
                    sys.stdout.write("".join(out_line) + "\n");
                    sys.stdout.flush();
            # For full updates, print the full line to the screen
            # For others, delete the "In progress..." column and update the same status line
            
            printWrite(globs['logfilename'], out_v, "".join(file_line));
            # Write the full line to the file.
        # The final status entry
        #####

    return cur_time;

#############################################################################

def errorOut(errnum, errmsg, globs):
# Formatting for error messages.
    fullmsg = "**Error " + str(errnum) + ": " + errmsg;
    border = "-" * len(fullmsg);
    fullstr = "\n" + border + "\n" + fullmsg + "\n" + border + "\n"
    printWrite(globs['logfilename'], globs['log-v'], "\n" + border + "\n" + fullmsg + "\n" + border + "\n");
    # Format and print the error message

    if globs['warnings'] != 0:
        warnmsg = "**Additionally there were " + str(globs['warnings']) + " warnings. Check the log file for more info";
        warnborder = border = "-" * len(warnmsg);
        printWrite(globs['logfilename'], globs['log-v'], "\n" + warnborder + "\n" + warnmsg + "\n" + warnborder + "\n");
    # Format and print information about warnings if there are any

    if globs['endprog']:
        globs['exit-code'] = 1;
        endProg(globs);
    else:
        printWrite(globs['logfilename'], globs['log-v'], "\nScript call: " + " ".join(sys.argv));
        sys.exit(1);
    # Exit the program

#############################################################################

def endProg(globs):
# A nice way to end the program.
    if globs['log-v'] == 1 :
        globs['log-v'] = 2;
    endtime = timeit.default_timer();
    totaltime = endtime - globs['starttime'];

    printWrite(globs['logfilename'], globs['log-v'], "# " + "=" * 175);
    printWrite(globs['logfilename'], globs['log-v'], "#\n# Done!");
    printWrite(globs['logfilename'], globs['log-v'], "# The date and time at the end is: " + getDateTime());
    printWrite(globs['logfilename'], globs['log-v'], "# Total execution time:            " + str(round(totaltime,3)) + " seconds.");
    printWrite(globs['logfilename'], globs['log-v'], "# Output directory for this run:   " + globs['outdir']);
    printWrite(globs['logfilename'], globs['log-v'], "# Log file for this run:           " + globs['logfilename']);

    if globs['warnings'] != 0:
        printWrite(globs['logfilename'], globs['log-v'], "\n# GRAMPA finished with " + str(globs['warnings']) + " WARNINGS -- check log file for more info");

    if globs['exit-code'] != 0:
        printWrite(globs['logfilename'], globs['log-v'], "#\n# ERROR: NON-ZERO EXIT STATUS.");
        printWrite(globs['logfilename'], globs['log-v'], "# ERROR: GRAMPA FINISHED WITH ERRORS.");
        printWrite(globs['logfilename'], globs['log-v'], "# ERROR: PLEASE CHECK THE LOG FILE FOR MORE INFO: " + globs['logfilename'] + "\n#");
    elif not globs['mul-opt'] and not globs['check-nums']:
        printWrite(globs['logfilename'], globs['log-v'], "# ----------------------------------------");
        if not globs['mul-input-flag']:
            if globs['min-num'] != 0:
                import grampa_lib.mul_tree as MT
                printWrite(globs['logfilename'], globs['log-v'], "# The MUL-tree with the minimum parsimony score is MT-" + str(globs['min-num']) + ":\t" \
                    + MT.mulPrint(globs['min-tree'][0], globs['min-tree'][2]));
            else:
                printWrite(globs['logfilename'], globs['log-v'], "# The tree with the minimum parsimony score is the singly-labled tree (ST):\t" + globs['min-tree'][0]);
            printWrite(globs['logfilename'], globs['log-v'], "# Score = " + str(globs['min-score']));
            printWrite(globs['logfilename'], globs['log-v'], "# ----------------------------------------");        

    #print("# " + "=" * 125);
    printWrite(globs['logfilename'], globs['log-v'], "# " + "=" * 175);
    printWrite(globs['logfilename'], globs['log-v'], "#");
    sys.exit(globs['exit-code']);

#############################################################################

def osCheck(test_cmd):
# For the tests script. Need to know if we are on Windows in order to pass the command correctly.
    import platform
    if platform.system == 'Windows':
        test_cmd = " ".join(test_cmd);
    return test_cmd;

#############################################################################

def testPrep():
# Prepares the test command and calls the tests script.
    t_path = os.path.join(os.path.dirname(__file__), "tests.py");
    pyver = sys.version[:3];
    try:
        python_cmd = "python" + pyver
        test_cmd = [python_cmd, t_path, python_cmd];
        subprocess.call(osCheck(test_cmd));
    except OSError:
        python_cmd = "python"
        test_cmd = [python_cmd, t_path, python_cmd];
        subprocess.call(osCheck(test_cmd));


#############################################################################




















































































































































































































































































































































































































































































































































































































































































































def simpson():
    s = """
                      @                                                           
                 CC   CQ                                                      
                /CCB @CC                                                      
            GCCS CCCCCCC7                                                     
             @CCCCCCCC@@@                                                     
            @@@CCCCCCCCCCCCC                                                  
            @CCCCCCCCCCCCCCCCC/                                               
              OCCCCCCCCCCCCCCCCC                                              
              @CCCCCCCCCCCCCCCCCCC@                                           
              CCCCCCCCCCC@QCCCCCC@CCC(                                        
              6CCCCCCCCCCCOCCCCCCCCC@es@                                      
              @CCCCCCCCCCCCCCCCCCCCCCCCCC                                     
              ^CCCCCCCCCCCCC@      @K      R                                  
               CCCC@CCCCCCC                                                   
              ~CCCC@CCCCCC@          G     //                   SCC@  @CC~    
              @CCCCsCCCCCC#    S@    #       RS@                CCCB 6CCC     
              @CCCCCCCCCOOG        S/@OC@CCSR  @/K             GCCCC@CCCC @CC 
               @CC@GCCCCCCCCS      K67@CCCCCCCCG @             CCCCKCCCC@CCCC 
              sCCCCCCCCCCCCCCCCCC77777SCCCCCCCCCCC        @@   CCCCKCCCCCCCCK 
              CCCGRCCCCCCCCCCCCCCC777@CC@CCCCCCCCC       @CCCC@CCCCCCCCCCCCR  
              @CCCK@CCCCCCCCCCCCCCCSQ(((((((((((@         6CCCCCCCCCCCCCCCC(  
               ^QCCCCCCCCCCCCCCCC%(((((((((((((((((%@s     @CCCCCCCCCCCCCCC   
                3CBCCCCCCCCCCR(((((((((((((((((((((((((%    CCCCCCC@CCCCCCC   
               /CCCCCCCCCCC(((((((((((((((((((@((((@((@     SCCCCCCCCCCCCCC   
               6CCCC@@CCCC(((((@   #@@@@       @@K          @KCCCCCC@CCCCCC   
               CCCCCBCCCCC(((((@@@@@@@@(((3                @((#CCCCCKCCCCC    
               CCCCSCCCCCC@((((@KK@KR@(((               %~(%(((%sCCCCCCCC     
               ~@CCCCCCC@CCC((((((O@@@%((C             @~(((@(((((~(6((%      
              (((((@CCCCRCCC@(((((((((((@@            @(((((((@t(((((((       
              ((((((((@@@@CCC#(((((((((@%@(@         s(((((((((((((Ct(        
             %(@(((((((((((((((((((@   /((((       @(((((((((((((((((@        
            sR(((@((((((((((((((((@      (#%@G((((((((((((((((((((((@         
          7~((((((((t@((((((((((G(s      ((((((((((((((((((((((((((R          
         B((((((((((((((s#((((((((((/  s(e((((((((((((((((((((((((R           
        @((((((((((((((((((((((((((((@((((@((((((((((((((((((((((/            
       /(((((((((((((((((((((((((((((((((((((((((((((((((((((((C              
       ~(((((((((((((((((((((((((((((e((((((((((((((((((((((((@               
      @((((((((((((((((((((((((((((((@(((@(#((((((((((((((((K                 
      G((((((((((((((((((((((((((((((s((((%@(((((((((((((((6                  
      (((((((((((((((((((((((((((((((((((((e(((((((((((Q/                     
     G((((((((((((((@((((((((((((((((((((((%(((((((@                          
     @(((((((((((((C(((((((((((((((((((t(((G(                                 
     ~(((((((((((((((((((((((((((((((((t(((e(s                                
    @(((((((((((((((((((((((((((((((((((((((((                                
    7(((((((((((((6((((((((((((((((((((Q((((((%   
    """

    print(s);
