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
        
    def getUpperOrLowerLimit(self, values, multiplier, limit_type):
        """
        """
        if limit_type not in ("upper", "lower"):
            raise ValueError(limit_type)
        mean = numpy.mean(values)
        std = numpy.std(values, ddof=1)
        if limit_type == "upper":
            return mean + (multiplier * std)
        elif limit_type == "lower":
            return mean - (multiplier * std)
        
    def plotSTDLimits(self, ax, values, vmin, vmax, no_std, line_color, plot_axis, plot_upper, plot_lower):
        """
        Plots three standard deviation lines to either X or Y axis. Both limits (upper and lower) or only one of
        those will be plotted.
        """
        if plot_axis not in ["x", "y"]:
            raise ValueError(plot_axis)
        upper = self.getUpperOrLowerLimit(values, no_std, "upper")
        lower = self.getUpperOrLowerLimit(values, no_std, "lower")
        if plot_upper:
            if plot_axis == "x":
                ax.plot([vmin, vmax], [upper,upper], line_color)
            else:
                ax.plot([upper, upper], [vmin, vmax], line_color)
        if plot_lower:
            if plot_axis == "x":
                ax.plot([vmin, vmax], [lower, lower], line_color)
            else:
                ax.plot([lower, lower], [vmin, vmax], line_color)
        return ax

        
    def plotThreeSTDLimits(self, ax, vals, vmin, vmax, plot_axis = "x", plot_upper = True, plot_lower = True):
        """
        Maps three standard deviation lines for given axis
        """
        stds_and_colors = ((5, "g--"), (4, "y--"), (3, "r--"))
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
        #color_map = ['green'] * len(x_vals)
        
        xmin = min(x_vals)
        xmax = max(x_vals)
        ymin = min(y_vals)
        ymax = max(y_vals)
        
        # plot the three standard deviation lines for each axis they are requested
        if plot_std_x_axis:
            ax = self.plotThreeSTDLimits(ax, y_vals, xmin, xmax, plot_axis = "x")
        if plot_std_y_axis:
            ax = self.plotThreeSTDLimits(ax, x_vals, ymin, ymax, plot_axis = "y", plot_lower = False)
        
        # calculate and plot the selection line for x axis (unit in std:s)
        if x_std_limit != None:
            x_limit = self.getUpperOrLowerLimit(x_vals, x_std_limit, "upper")
            ax.plot([x_limit, x_limit], [ymin, ymax], 'b--')
            #color_map = self.getColorMap(x_vals, x_std_limit, color_map, limit_type = "upper", use_std_limit = True)
        
        # calculate and plot the selection line for y axis lower limit (unit in std:s)
        if y_std_limit != None:
            y_limit = self.getUpperOrLowerLimit(y_vals, y_std_limit, "lower")
            ax.plot([xmin, xmax], [y_limit, y_limit], 'b--')
            #color_map = self.getColorMap(y_vals, y_std_limit, color_map, use_std_limit = True)
            
        # calculate and plot the selection line for y axis upper limit (unit in std:s)
        if y_std_limit != None:
            y_limit = self.getUpperOrLowerLimit(y_vals, y_std_limit, "upper")
            ax.plot([xmin, xmax], [y_limit, y_limit], 'b--')
            #color_map = self.getColorMap(y_vals, y_std_limit, color_map, limit_type = "upper", use_std_limit = True)
            
        # plot x axis limit
        if x_val_limit != None:
            ax.plot([x_val_limit, x_val_limit], [min(y_vals), max(y_vals)], 'b--')
            #color_map = self.getColorMap(x_vals, x_val_limit, color_map, limit_type = "upper")
            
        # plot y axis lower limit 
        if y_val_limit != None:
            ax.plot([min(x_vals), max(x_vals)], [y_val_limit, y_val_limit], 'b--')
            #color_map = self.getColorMap(y_vals, y_val_limit, color_map)
        
        # plot y axis upper limit 
        if y_val_limit != None:
            ax.plot([min(x_vals), max(x_vals)], [y_val_limit, y_val_limit], 'b--')
            #color_map = self.getColorMap(y_vals, y_val_limit, color_map, limit_type = "upper")
            
        # use log scale if requested
        if x_log_scale:
            ax.set_xscale('log')
        if y_log_scale:
            ax.set_yscale('log')
        #ax.scatter(x_vals, y_vals, c=color_map)
        ax.scatter(x_vals, y_vals, c=col_map)
        #print ("\n***********************************************")
        #print ("Plotter: given col map red: %d plotter red: %d" % (col_map.count("red"), color_map.count("red")))
        canvas.print_figure(plot_fn)
        #return color_map
        return
    
    
    def getColorMap(self, values, limit, col_map, limit_type = "lower", use_std_limit = False):
        """
        Makes color map based on given limits. If standard deviation based limit is given, new limit is calculated,
        otherwise the value given is used directly. Any values smaller / larger than the limit is marked red.
        """
        if limit_type not in ("lower", "upper"):
            raise ValueError(limit_type)
        # calculate limit if std dev is used
        if use_std_limit:
            if limit_type == "upper":
                limit = self.getUpperOrLowerLimit(values, limit, "upper")
            else: 
                limit = self.getUpperOrLowerLimit(values, limit, "lower")
        for i, val in enumerate(values):
            if limit_type == "upper":
                if val > limit:
                    col_map[i] = "red"
            elif limit_type == "lower":
                if val < limit:
                    col_map[i] = "red"
            else:
                raise ValueError(limit_type)
        return col_map

        
#     def plotHistogram(self, values, no_bins = 50):
#         print ("Printing histogram of %d values" % len(values))
#         values = numpy.array(values)
#         if numpy.isnan(values).any():
#             print ("Found numpy nans!")
#             values = values[~numpy.isnan(values)]
#             #import math
#             #values = [value for value in values if not math.isnan(value)]
#             print ("No entries in cleaned list: %s" % len(values))
#         plt.hist(values, no_bins)
#         plt.show()
        