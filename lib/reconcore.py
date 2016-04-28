#############################################################################
# Reconciliation CORE functions
# Gregg Thomas
# August 2013-present
# Forked from core on 12.13.2015
#############################################################################

import sys, datetime

#############################################################################

def errorOut(errnum, errmsg):
# Formatting for error messages.
	fullmsg = "|**Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	print "\n" + border + "\n" + fullmsg + "\n" + border + "\n";

#############################################################################

def getLogTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y-%I.%M.%S");

#############################################################################

def getDateTime():
# Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

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
	if v == 1 or v == -2:
		print outline;
	f = open(o_name, "a");
	f.write(outline + "\n");
	f.close();

#############################################################################

def hErrorCheck(spec_check, sinfo, spec_type, hybrid_clades, copy_clades):

	if hybrid_clades and not all(h in sinfo for hybrid_list in hybrid_clades for h in hybrid_list if not h.isdigit()):
		errorOut(8, "Not all hybrid species (-h1) are present in your species tree!");
		return 1

	if hybrid_clades and not all("<" + h + ">" in sinfo for hybrid_list in hybrid_clades for h in hybrid_list if h.isdigit()):
		errorOut(9, "Not all hybrid nodes (-h1) are present in your species tree!");
		return 1

	if copy_clades and not all(c in sinfo for copy_list in copy_clades for c in copy_list if not c.isdigit()):
		errorOut(10, "Not all copy species (-h2) are present in your species tree!");
		return 1;

	if copy_clades and not all("<" + c + ">" in sinfo for copy_list in copy_clades for c in copy_list if c.isdigit()):
		errorOut(11, "Not all copy nodes (-h2) are present in your species tree!");
		return 1;

	if spec_type == 's' and any(spec_check.count(n) > 1 for n in spec_check):
		errorOut(12, "You have entered a tree type (-t) of 's' but there are labels in your tree that appear more than once!");
		return 1;

	if spec_type == 'm' and any(spec_check.count(h) not in [1,2] for h in spec_check):
		errorOut(13, "You have entered a tree type (-t) of 'm', species in your tree should appear exactly once or twice.");
		return 1;

	return 0;

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

	print s;




