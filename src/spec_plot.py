from bokeh.plotting import figure
from src.mycolors import mycolors

def spec_plot(datalist,colorlist=None,**kwargs):
    golden=1.61803
    height=400
    #plot_width=int(height*golden)
    if 'plot_height' not in kwargs:
        kwargs['plot_height']= height
        kwargs['plot_width']=int(height*golden)
    if 'plot_width' not in kwargs:
        kwargs['plot_width']=int(kwargs['plot_height']*golden)
        
    myplot=figure(**kwargs)
    myplot.outline_line_alpha=0
    
    #Xaxis prefs
    myplot.xgrid.grid_line_color = None
    myplot.xaxis.axis_label="Evolution Time (μs)"
    myplot.xaxis.axis_label_text_font='helvetica'
    myplot.xaxis.axis_label_text_font_style='normal'
    myplot.xaxis.axis_label_text_font_size='24pt'
    myplot.xaxis.major_label_text_font='helvetica'
    myplot.xaxis.major_label_text_font_size='16pt'

    #Yaxis prefs
    myplot.ygrid.grid_line_alpha = 0.5
    myplot.ygrid.grid_line_dash = [6, 4]
    myplot.yaxis.major_label_text_font='helvetica'
    myplot.yaxis.major_label_text_font_size='16pt'

    if not colorlist: #if colorlist is not present, use default scheme
        for i in range(int(len(datalist)/2)):
            myplot.line(datalist[i*2],datalist[i*2+1], color=mycolors(i), line_width=1)
        return myplot
    else:
        for i in range(int(len(datalist)/2)):
            myplot.line(datalist[i*2],datalist[i*2+1], color=colorlist[i], line_width=1)
        return myplot