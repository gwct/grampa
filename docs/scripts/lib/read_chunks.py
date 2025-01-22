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

def readNav(active_url, path1, path2, path3):
    navfile = "../html-chunks/nav.html";
    navlines = open(navfile, "r").readlines()
    for x in range(len(navlines)):
        if active_url in navlines[x]:
            navlines[x] = navlines[x].replace(active_url, "#");
            if 'class="nav-link"' in navlines[x]:
                navlines[x] = navlines[x].replace('class="nav-link"', 'class="nav-link" id="active"');
        if "PATH_STR_REPLACE1" in navlines[x]:
            navlines[x] = navlines[x].replace("PATH_STR_REPLACE1", path1);
        if "PATH_STR_REPLACE2" in navlines[x]:
            navlines[x] = navlines[x].replace("PATH_STR_REPLACE2", path2);
        if "PATH_STR_REPLACE3" in navlines[x]:
            navlines[x] = navlines[x].replace("PATH_STR_REPLACE3", path3);         
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
    year = datetime.now().strftime("%Y");
    zone = subprocess.check_output("date +%Z").decode().strip();
    return open(footerfile, "r").read().replace("CURYEAR", year).replace("DATETIME", now + " " + zone);