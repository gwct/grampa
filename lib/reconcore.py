#############################################################################
# Reconciliation CORE functions
# Gregg Thomas
# August 2013-present
# Forked from core on 12.13.2015
#############################################################################

import sys, os, subprocess, datetime, time, opt_parse as OP

#############################################################################

def errorOut(errnum, errmsg):
# Formatting for error messages.
	OP.optParse(1);
	fullmsg = "|**Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	print("\n" + border + "\n" + fullmsg + "\n" + border + "\n");
	sys.exit();

#############################################################################

def endProg(starttime, outfilename, main_v):
		endtime = time.time();
		totaltime = endtime - starttime;
		printWrite(outfilename, main_v, "# LOG: Total execution time: " + str(round(totaltime,3)) + " seconds.");
		printWrite(outfilename, main_v, "# =========================================================================");
		sys.exit();

#############################################################################

def getLogTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y-%I.%M.%S");

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

#############################################################################

def filePrep(filename, header=""):
# Function to initialize output files with headers or as blank files.
	if header != "" and header[-1] != "\n":
		header += "\n";
	outfile = open(filename, "w");
	outfile.write(header);
	outfile.close();

#############################################################################

def loadingBar(counter, length, done, bars):
#This function serves as a text loading bar for long scripts with counters. The following
#lines must be added within the script to initialize and terminate the script:
#Initilization:
#numlines = core.getFileLen(alnfilename);
#numbars = 0;
#donepercent = [];
#i = 0;
#Termination:
#	pstring = "100.0% complete.";
#	sys.stderr.write('\b' * len(pstring) + pstring);
#	print "\nDone!";
#
#If length is lines in a file use the core.getFileLen function to get the number of lines in the file

	percent = float(counter) / float(length) * 100.0;
	percentdone = int(percent);

	p = str(percent)
	pstring = " " + p[:5] + "% complete.";

	if percentdone % 2 == 0 and done != None and percentdone not in done:
		loading = "";
		loading = "[";
		j = 0;
		while j <= bars:
			loading = loading + "*";
			j = j + 1;
		while j < 50:
			loading = loading + "-";
			j = j + 1;
		loading = loading + "]";

		loading = loading + "                 ";
		sys.stderr.write('\b' * len(loading) + loading);

		done.append(percentdone);
		bars = bars + 1;

	sys.stderr.write('\b' * len(pstring) + pstring);

	return bars, done;

#############################################################################

def printWrite(o_name, v, o_line1, o_line2="", pad=0):
#Function to print a string AND write it to the file.
	if o_line2 == "":
		outline = o_line1;
	else:
		outline = o_line1 + " "*(pad-len(o_line1)) + o_line2;
	if v in [1,2]:
		print(outline);
	f = open(o_name, "a");
	f.write(outline + "\n");
	f.close();

#############################################################################

def printStep(step, msg, v):
	if v != -1:
		print msg;
	return step+1;

#############################################################################
def osCheck(test_cmd):
	import platform
	if platform.system == 'Windows':
		test_cmd = " ".join(test_cmd);
	print test_cmd;
	return test_cmd;

#############################################################################

def testPrep():
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

## Old code for the ST branch output.

# From reconLCA:
	# node_counts = {};
	# for node in sinfo:
	# 	node_counts[node] = [0,0];

	# node_counts[lca_maps[d1][0]][1] += d1_loss;
	# node_counts[lca_maps[d2][0]][1] += d2_loss;

	# for node in lca_maps:
	# 	if dups[node] != 0:
	# 		node_counts[lca_maps[node][0]][0] += 1;
	# Place the duplications on their maps in the species tree.

# From stdLCA:
	# branch_outline = "\t";
	# for node in sorted(st_node_counts.keys()):
	# 	tot_node_counts[node][0] += st_node_counts[node][0];
	# 	tot_node_counts[node][1] += st_node_counts[node][1];

	# 	branch_outline += "\t" + node + ":" + str(st_node_counts[node][0]) + "," + str(st_node_counts[node][1])
	# RC.printWrite(detoutfilename, v, branch_outline);
	# Write the branch gain/loss scores.
# branch_outline = "Total branch scores for ST:\t" + "\t".join([node + ":" + str(tot_node_counts[node][0]) + "," + str(tot_node_counts[node][1]) for node in sorted(tot_node_counts.keys())]);
# RC.printWrite(detoutfilename, v, branch_outline);

# From grampa.py:
	# branch_outline = "\t";
	# for node in sorted(mt_node_counts.keys()):
	# 	tot_node_counts[node][0] += mt_node_counts[node][0];
	# 	tot_node_counts[node][1] += mt_node_counts[node][1];

	# 	branch_outline += "\t" + node + ":" + str(mt_node_counts[node][0]) + "," + str(mt_node_counts[node][1])
	# RC.printWrite(detoutfilename, v, branch_outline);
	# Write the branch gain/loss scores.

# branch_outline = "Total branch scores for MT" + str(mul_num) + ":\t" + "\t".join([node + ":" + str(tot_node_counts[node][0]) + "," + str(tot_node_counts[node][1]) for node in sorted(tot_node_counts.keys())]);
# RC.printWrite(detoutfilename, v, branch_outline);

#############################################################################
## The old loss counting function.

# def mulLossCount(lc_ginfo, lc_minfo, lc_maps, lc_dups, lc_node_counts):
# # Given two trees (dictionaries), a mapping between them, and the duplication nodes,
# # this function counts the number of losses. Depths are 0 based.

# 	loss_count = 0;
# 	for g in lc_ginfo:
# 		if lc_ginfo[g][2] == 'root':
# 			glosses = len(RT.nodeDepth(lc_maps[g][0],lc_minfo));
# 		# The number of losses at the root of the gene tree is equal to the depth of its map.

# 		else:
# 			curanc = lc_ginfo[g][1];
# 			ancdepth = len(RT.nodeDepth(lc_maps[curanc][0],lc_minfo));
# 			gdepth = len(RT.nodeDepth(lc_maps[g][0],lc_minfo));

# 			glosses = 0;
# 			glosses = gdepth - ancdepth - 1;

# 			if lc_dups[curanc] != 0:
# 				glosses = glosses + 1;

# 		if glosses != 0:
# 			loss_count += glosses;
# 			lc_node_counts[lc_maps[g][0]][1] += glosses; 	

# 	# for m in lc_minfo:
# 	# 	print m;
# 	# 	if lc_minfo[m][2] == 'root' and [m] not in list(lc_maps.values()):
# 	# 		loss_count = loss_count + 1;
# 	# 		break;
# 	# Accounts for cases where h2 puts one clade at the root of the MUL-tree
# 	#print loss_count;
# 	#sys.exit();
# 	return loss_count, lc_node_counts;

#############################################################################













































































































































































































































































































































































































































































































































































































































































































































































































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

#############################################################################