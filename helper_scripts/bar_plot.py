from plotly.offline import plot
import plotly.graph_objs as go
import plotly.plotly as py
# Option parsing and import of plot libraries if no errors.

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

def sideBarPlot(traces,maintitle,outname,barcol='rgb(0,102,51)',plotcol='#e1e1ea',bgcol='#fffae6',w=1000,h=1000,bmar=150):
	#data = [go.Bar(x=xdata,y=ydata,marker=dict(color=barcol),opacity=0.6)];
	data = [];
	for trace in traces:
		cur_trace = go.Bar(x=trace[0],y=trace[1],name=trace[2]);
		data.append(cur_trace);


	layout = go.Layout(
		barmode='group',
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
