import sys, os

############################################

def barPlot(xdata,ydata,xtitle,ytitle,maintitle,outname,barcol='rgb(0,102,51)',plotcol='#e1e1ea',bgcol='#fffae6',w=1000,h=1000,bmar=150):
	data = [go.Bar(x=xdata,y=ydata,marker=dict(color=barcol),opacity=0.6)];
	layout = go.Layout(
		autosize=False,
		width=w,
		height=h,
		paper_bgcolor=bgcol,
		plot_bgcolor=plotcol,
		title=maintitle,
		titlefont=dict(
			family="Arial, sans-serif",
			size=30,
			),
		xaxis=dict(
			title=xtitle,
			titlefont=dict(
				family="Arial, sans-serif",
				size=20,
				color="#737373"
				),
			),
		yaxis=dict(
			title=ytitle,
			titlefont=dict(
				family="Arial, sans-serif",
				size=20,
				color="#737373"
				)
			)
		);
	fig = go.Figure(data=data, layout=layout);
	plot(fig, filename=outname);

############################################

def scatterPlot(xdata,ydata,xtitle,ytitle,maintitle,outname,barcol='rgb(0,102,51)',plotcol='#e1e1ea',bgcol='#fffae6',w=1000,h=500,bmar=150):
	data = [go.Scatter(x=xdata,y=ydata,mode='markers',opacity=0.6)];
	layout = go.Layout(
		autosize=False,
		width=w,
		height=h,
		margin=go.Margin(
        	l=70,
        	r=20,
        	b=150,
        	t=70,
        	pad=0
    	),
		paper_bgcolor=bgcol,
		plot_bgcolor=plotcol,
		title=maintitle,
		titlefont=dict(
			family="Arial, sans-serif",
			size=30,
			),
		xaxis=dict(
			title=xtitle,
			titlefont=dict(
				family="Arial, sans-serif",
				size=20,
				color="#737373",
				),
			tickangle=90
			),
		yaxis=dict(
			title=ytitle,
			titlefont=dict(
				family="Arial, sans-serif",
				size=20,
				color="#737373"
				)
			)
		);
	fig = go.Figure(data=data, layout=layout);
	plot(fig, filename=outname);

############################################

if len(sys.argv) != 3 or "-h" in sys.argv:
	print("\n# Usage: grampa_plot.py [input file] [output file]");
	print("# ---> [input file] must be a grampa output file.")
	print("# ---> [output file] will be an html file with your plot.\n")
	sys.exit();

infilename = sys.argv[1];
outfilename = sys.argv[2];
if outfilename[len(outfilename)-5:] != ".html":
	outfilename += ".html";

from plotly.offline import plot
import plotly.graph_objs as go
import plotly.plotly as py
# Option parsing and import of plot libraries if no errors.

score_dict = {};

for line in open(infilename):
	if line[0] == "#":
		continue;

	line = line.strip().split("\t");
	score_dict[line[1] + "-" + line[2]] = int(line[4]);

sorted_keys = sorted(score_dict, key=score_dict.get)
sorted_vals = [];

max_len = -999;
for key in sorted_keys:
	sorted_vals.append(score_dict[key]);

	if len(key) > max_len:
		max_len = len(key);

bot_margin = max_len * 15;

scatterPlot(sorted_keys,sorted_vals,"H1-H2 Node", "Score", "GRAMPA Results: " + os.path.basename(infilename), outfilename, bmar=bot_margin);
