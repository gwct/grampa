############################################################
# For murine web development, 11.17
# Functions to read static html chunks
############################################################

def readHead(title, active_url):
    page_name = active_url.replace(".html", "");
    headfile = "../html-chunks/head.html";
    return open(headfile, "r").read().replace("TMPTITLE", title).replace("PAGECSS", page_name);

def readResultsHead(title):
    headfile = "../html-chunks/results_head.html";
    return open(headfile, "r").read().replace("TMPTITLE", title);

def readNav(active_url):
    navfile = "../html-chunks/nav.html";
    navlines = open(navfile, "r").readlines()
    for x in range(len(navlines)):
        if active_url in navlines[x]:
            navlines[x] = navlines[x].replace(active_url, "#");
            if 'class="nav_link"' in navlines[x]:
                navlines[x] = navlines[x].replace('class="nav_link"', 'class="nav_link" id="active"');
    return "".join(navlines);

def readWheatNav(runtype):
    navfile = "../html-chunks/wheat_nav.html";
    return open(navfile, "r").read().replace("TMPTYPE", runtype);

def readYeastNav():
    navfile = "../html-chunks/yeast_nav.html";
    return open(navfile, "r").read();

def readFooter():
    import time, subprocess
    from datetime import datetime
    footerfile = "../html-chunks/footer.html";
    now = datetime.now().strftime("%m/%d/%Y %H:%M:%S");
    zone = subprocess.check_output("date +%Z").decode().strip();
    return open(footerfile, "r").read().replace("DATETIME", now + " " + zone);