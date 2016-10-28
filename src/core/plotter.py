'''
Created on 22.12.2015

@author: paakkone
'''
import os
import numpy
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

class Plotter(object):
    """
    Does the plotting for the analysis
    """
    def __init__(self):
        pass
        
    def getUpperAndLowerLimit(self, values, multiplier):
        """
        Returns  upper and lower limits of given range (average +- std * std_multiplier)
        """
        mean = numpy.mean(values)
        std = numpy.std(values, ddof=1)
        upper = mean + (multiplier * std)
        lower = mean - (multiplier * std)
        return lower, upper
       
        
    def plotSTDLimits(self, ax, values, vmin, vmax, no_std, line_color, plot_axis, plot_upper, plot_lower):
        """
        Plots three standard deviation lines to either X or Y axis. Both limits (upper and lower) or only one of
        those will be plotted.
        """
        if plot_axis not in ["x", "y"]:
            raise ValueError(plot_axis)
        lower, upper = self.getUpperAndLowerLimit(values, no_std)
        if plot_upper:
            if plot_axis == "x":
                ax.axhline(y=upper, color = line_color, linestyle = "dashed")
            else:
                ax.axvline(x=upper, color=line_color, linestyle = "dashed")
        if plot_lower:
            if plot_axis == "x":
                ax.axhline(y=lower, color = line_color, linestyle = "dashed")
            else:
                ax.axvline(x=lower, color=line_color, linestyle = "dashed")
        return ax

        
    def plotThreeSTDLimits(self, ax, vals, vmin, vmax, plot_axis = "x", plot_upper = True, plot_lower = True):
        """
        Maps three standard deviation lines for given axis
        """
        stds_and_colors = ((5, "g"), (4, "y"), (3, "r"))
        for std, col in stds_and_colors:
            ax = self.plotSTDLimits(ax, vals, vmin, vmax, std, col, plot_axis, plot_upper, plot_lower)
        return ax

    def plotScatterPlot(self, plot_fn, x_vals, y_vals, x_val_limit = None, y_val_limit = None, x_std_limit = None, 
                        y_std_limit = None, figure_title = "Title", x_title = "x-axis", 
                        y_title = "y_axis", x_log_scale = False, y_log_scale = False, 
                        plot_std_x_axis = True, plot_std_y_axis = False, col_map = None):
        """
        Makes scatter plot and saves it into png file
        """
        # setup the plot and its titles and axis labels
        fig = Figure()
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_subplot('111')
        ax.set_title(figure_title)
        ax.set_xlabel(x_title)
        ax.set_ylabel(y_title)
        
        xmin = min(x_vals)
        xmax = max(x_vals)
        ymin = min(y_vals)
        ymax = max(y_vals)
        
        # plot the three standard deviation lines for each axis they are requested
        if plot_std_x_axis:
            ax = self.plotThreeSTDLimits(ax, y_vals, xmin, xmax, plot_axis = "x")
        if plot_std_y_axis:
            ax = self.plotThreeSTDLimits(ax, x_vals, ymin, ymax, plot_axis = "y")
            #ax = self.plotThreeSTDLimits(ax, x_vals, ymin, ymax, plot_axis = "y", plot_lower = False)
        
        # calculate and plot the selection line for x axis (unit in std:s)
        if x_std_limit != None:
            x_lower, x_upper = self.getUpperAndLowerLimit(x_vals, x_std_limit)
            ax.axvline(x=x_lower, color='b')
            ax.axvline(x=x_upper, color='b')

        # calculate and plot the selection line for y axis lower limit (unit in std:s)
        if y_std_limit != None:
            y_lower, y_upper = self.getUpperAndLowerLimit(y_vals, y_std_limit)
            ax.axhline(y=y_lower, color = 'b')
            ax.axhline(y=y_upper, color = 'b')
            
        # plot x axis limit
        if x_val_limit != None:
            ax.axvline(x=x_val_limit, color='b')
            
        # plot y axis lower limit 
        if y_val_limit != None:
            ax.axhline(y=y_val_limit, color = 'b')
            
        # use log scale if requested
        if x_log_scale:
            ax.set_xscale('log')
        if y_log_scale:
            ax.set_yscale('log')
        ax.scatter(x_vals, y_vals, c=col_map)
        canvas.print_figure(plot_fn)
        return
    