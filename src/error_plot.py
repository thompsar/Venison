from bokeh.palettes import Colorblind6 as palette
from bokeh.plotting import figure,show
from bokeh.layouts import column, row, gridplot
import numpy as np

def error_plot(landscape,statistics,cutoffs = [0.95,0.75,0.6], refresh = False):

    __, npop, nparams = landscape.shape
    
    error_surface_plot = figure(plot_width = 600, 
                                plot_height=600,
                                x_range = (2,7)
#                                 x_range = (
#                                     np.floor(landscape[:,:,0].min()),
#                                     np.ceil(landscape[:,:,0].max())
#                                 )
                               )

    error_surface_plot.xaxis.axis_label="Gaussian Center (nm)"
    error_surface_plot.xaxis.axis_label_text_font_size='14pt'
    error_surface_plot.yaxis.axis_label="Gaussian Width (nm)"
    error_surface_plot.yaxis.axis_label_text_font='helvetica'
    error_surface_plot.yaxis.axis_label_text_font_style='normal'
    error_surface_plot.yaxis.axis_label_text_font_size='14pt'
    error_surface_plot.yaxis.major_label_text_font_size='14pt'
    error_surface_plot.xaxis.major_label_text_font_size='14pt'
    
    error_surface_plot.ygrid.grid_line_alpha = 0.5
    error_surface_plot.xgrid.grid_line_alpha = 0.5


    for color, cutoff in zip([5,2,0],cutoffs):
        for pop in range(npop):
            error_surface_plot.circle(landscape[statistics <= cutoff ,pop,0],landscape[statistics <= cutoff,pop,1],color = palette[color],size=1)

        
    def slice_plot (xs, statistics, xlabel,ylabel,title,**kwargs):
    #     plot_width = plot_width, plot_height=plot_height
        sub_plot = figure(**kwargs)
        sub_plot.title.text = title
        sub_plot.xaxis.axis_label=xlabel
        sub_plot.xaxis.axis_label_text_font_size='10pt'
        sub_plot.xaxis.major_label_text_font_size='10pt'
        sub_plot.yaxis.major_label_text_font_size='10pt'
        sub_plot.yaxis.axis_label=ylabel
        sub_plot.yaxis.axis_label_text_font='helvetica'
        sub_plot.yaxis.axis_label_text_font_style='normal'
        sub_plot.yaxis.axis_label_text_font_size='10pt'
        sub_plot.ygrid.grid_line_alpha = 0.5
        sub_plot.xgrid.grid_line_alpha = 0.5
        sub_plot.xaxis.minor_tick_line_color = None
        sub_plot.yaxis.minor_tick_line_color = None
    #     left,rightkwargs['x_range']
        for color, cutoff in zip([5,2,0],[0.95,0.75,0.6]):
            sub_plot.circle(xs[statistics <= cutoff],statistics[statistics <= cutoff],color = palette[color],size=1)
        sub_plot.line(np.linspace(*kwargs['x_range']),[0.5]*50 ,color = 'black', line_dash = 'dashed')
        return sub_plot

    axes_labels = ['r (nm)','σ (nm)','Mole Fraction']
#     xranges = ((np.floor(landscape[:,:,0].min()),np.ceil(landscape[:,:,0].max())),(0,1),(0,1))
    xranges = xranges = ((2,7),(0,1),(0,1))
    yranges = (0,1)

    surface = gridplot([[error_surface_plot]])
    
    if npop>1:
    
        slicegrid = gridplot([[slice_plot(landscape[:,j,i],statistics,axes_labels[i],'F-ratio','Population '+str(j+1), plot_width=200, plot_height=200, x_range=xranges[i], y_range=yranges) for i in range(nparams)] for j in range(npop)]) 
    else:
        slicegrid = gridplot([[slice_plot(landscape[:,j,i],statistics,axes_labels[i],'F-ratio','Population '+str(j+1), plot_width=200, plot_height=200, x_range=xranges[i], y_range=yranges) for i in range(nparams-1)] for j in range(npop)]) 
    
    if not refresh:
    
        full_plot = show(row(surface,slicegrid),notebook_handle = True)

        return full_plot, surface, slicegrid
    
    return surface.children, slicegrid.children