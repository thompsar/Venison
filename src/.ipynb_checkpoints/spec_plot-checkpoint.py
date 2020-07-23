from bokeh.plotting import figure
from src.mycolors import mycolors
import numpy as np
from bokeh.palettes import Category10_10 as palette
import itertools  


##Need to reimplement colorlist for certain uses? (post implementing pallete)
def spec_plot(xs = None ,ys = None ,colorlist=None,**kwargs):
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

    colors = itertools.cycle(palette)  
    
    #this needs some serious cleaning/reworking
    if xs is not None and ys is not None:
    
        if (type(xs) == list) & (type(ys) == list):
            for i, color in zip(range(len(xs)),colors):
                myplot.line(xs[i],ys[i],line_width = 1, color = color)
        elif (type(xs) == np.ndarray) & (type(ys) == np.ndarray):
            try:
                for i, color in zip(range(xs.T.shape[1]),colors):
                    myplot.line(xs[i,:],ys[i,:],line_width = 1, color = color)
            except IndexError:
                    myplot.line(xs,ys,line_width = 1, color = next(colors))
    return myplot
    