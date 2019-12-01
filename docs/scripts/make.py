import sys, os, argparse

print()
print("###### Build site pages ######");
print("PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])))
print("# Script call: " + " ".join(sys.argv) + "\n----------");

parser = argparse.ArgumentParser(description="Gets stats from a bunch of abyss assemblies.");
parser.add_argument("--all", dest="all", help="Build all pages", action="store_true", default=False);
parser.add_argument("--index", dest="index", help="Without --all: build index.html. With --all: exlude index.html", action="store_true", default=False);
parser.add_argument("--about", dest="about", help="Without --all: build about.html. With --all: exlude about.html", action="store_true", default=False);
parser.add_argument("--readme", dest="readme", help="Without --all: build readme.html. With --all: exlude readme.html", action="store_true", default=False);
parser.add_argument("--example1", dest="example1", help="Without --all: build example1.html. With --all: exlude example1.html", action="store_true", default=False);
parser.add_argument("--example2", dest="example2", help="Without --all: build example2.html. With --all: exlude example2.html", action="store_true", default=False);
parser.add_argument("--example3", dest="example3", help="Without --all: build example3.html. With --all: exlude example3.html", action="store_true", default=False);
parser.add_argument("--yeast", dest="yeast", help="Without --all: build yeast.html. With --all: exlude yeast.html", action="store_true", default=False);
parser.add_argument("--yeastohno", dest="yeastohno", help="Without --all: build yeast_ohno.html. With --all: exlude yeast_ohno.html", action="store_true", default=False);
parser.add_argument("--wheatall", dest="wheatall", help="Without --all: build wheat_all.html. With --all: exlude wheat_all.html", action="store_true", default=False);
parser.add_argument("--wheatab", dest="wheatab", help="Without --all: build wheat_ab.html. With --all: exlude wheat_ab.html", action="store_true", default=False);
parser.add_argument("--wheatad", dest="wheatad", help="Without --all: build wheat_ad.html. With --all: exlude wheat_ad.html", action="store_true", default=False);
parser.add_argument("--wheatbd", dest="wheatbd", help="Without --all: build wheat_bd.html. With --all: exlude wheat_bd.html", action="store_true", default=False);
parser.add_argument("--performance", dest="performance", help="Without --all: build performance.html. With --all: exlude performance.html", action="store_true", default=False);
parser.add_argument("--links", dest="links", help="Without --all: build links.html. With --all: exlude links.html", action="store_true", default=False);
args = parser.parse_args();
# Input options.

#cwd = os.getcwd();
os.chdir("generators");

pages = {
    'index' : args.index,
    'about' : args.about,
    'readme' : args.readme,
    'example1' : args.example1,
    'example2' : args.example2,
    'example3' : args.example3,
    'yeast' : args.yeast,
    'yeastohno' : args.yeastohno,
    'wheatall' : args.wheatall,
    'wheatab' : args.wheatab,
    'wheatad' : args.wheatad,
    'wheatbd' : args.wheatbd,
    'performance' : args.performance,
    'links' : args.links,
}

if args.all:
    pages = { page : False if pages[page] == True else True for page in pages };

if pages['index']:
    os.system("python index_generator.py");

if pages['about']:
    os.system("python about_generator.py");

if pages['readme']:
    os.system("python readme_generator.py");

if pages['example1']:
    os.system("python example1_generator.py");

if pages['example2']:
    os.system("python example2_generator.py");

if pages['example3']:
    os.system("python example3_generator.py");

if pages['yeast']:
    os.system("python yeast_generator.py");

if pages['yeastohno']:
    os.system("python yeast_ohno_generator.py");

if pages['wheatall']:
    os.system("python wheat_all_generator.py");

if pages['wheatab']:
    os.system("python wheat_ab_generator.py");

if pages['wheatad']:
    os.system("python wheat_ad_generator.py");

if pages['wheatbd']:
    os.system("python wheat_bd_generator.py");

if pages['performance']:
    os.system("python performance_generator.py");

if pages['links']:
    os.system("python links_generator.py");
    
print("----------\nDone!");


