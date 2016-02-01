'''
Created on 22.12.2015

@author: paakkone
'''
import os
import numpy
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt

#print ("Matplotlib backend: %s" % matplotlib.get_backend())

class Plotter(object):
    """
    Does the plotting for the analysis
    """
    def __init__(self):
        pass
    
            
#     def getMADLimits(self, x, no_stds = 3, mad_to_std = False, double = False):
#         """
#         Calculates limits using median-absolute-deviation
#         """
#         x = numpy.array(x)
#         # check that less than 50 % of values are median
#         unique, counts = numpy.unique(x, return_counts = True)
#         max_count = max(counts)
#         #print ("No values: %d,  max_no_same_value: %d" % (len(x), max_count))
#         if max_count >= (0.5 * len(x)):
#             raise ValueError
#         median = numpy.median(x)
#         if double:
#             deviations = x - median
#             upper_devs = [dev for dev in deviations if dev >= 0]
#             upper_mad = numpy.median(upper_devs)
#             lower_devs = [dev for dev in deviations if dev <= 0]
#             lower_mad = numpy.median(lower_devs)
#             #print ("median: %s, no_x: %s, no_abs_dev: %s" % (median, len(x), len(abs_dev)))
#             upper = median + (no_stds * upper_mad)
#             lower = median + (no_stds * lower_mad)
#             print ("Upper MAD: %s, lower MAD: %s" % (upper_mad, lower_mad))
#         else:
#             abs_dev = abs(x - median)
#             if mad_to_std:
#                 mad = 1.4826 * numpy.median(abs_dev)
#             else:
#                 mad = numpy.median(abs_dev)
#             upper = median + (no_stds * mad)
#             lower = median - (no_stds * mad)
#             print ("MAD: %s" % mad)
#         print ("Median: %s, MAD upper: %s, lower %s" % (median, upper, lower))
#         return upper, lower
        
    def plotSTDLimits(self, ax, values, vmin, vmax, no_std, line_color, plot_axis, plot_upper, plot_lower):
        """
        Plots three standard deviation lines to either X or Y axis. Both limits (upper and lower) or only one of
        those will be plotted.
        """
        if plot_axis not in ["x", "y"]:
            raise ValueError(plot_axis)
        mean = numpy.mean(values)
        std = numpy.std(values, ddof=1)
        upper = mean + (no_std * std)
        lower = mean - (no_std * std)
        if plot_upper:
            if plot_axis == "x":
                ax.plot([vmin, vmax], [upper,upper], line_color)
                #print ("plotting [%3.3f - %3.3f], [%3.3f - %3.3f] " % (vmin, vmax, upper, upper))
            else:
                ax.plot([upper, upper], [vmin, vmax], line_color)
                #print ("plotting [%3.3f - %3.3f], [%3.3f - %3.3f] " % (upper, upper, vmin, vmax))
        if plot_lower:
            if plot_axis == "x":
                ax.plot([vmin, vmax], [lower, lower], line_color)
                #print ("plotting [%3.3f - %3.3f], [%3.3f - %3.3f] " % (vmin, vmax, lower, lower))
            else:
                ax.plot([lower, lower], [vmin, vmax], line_color)
                #print ("plotting [%3.3f - %3.3f], [%3.3f - %3.3f] " % (lower, lower, vmin, vmax))
#             ax.plot([min(x_values), max(x_values)], [lower,lower], line_color)
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
                        plot_std_x_axis = True, plot_std_y_axis = False):
        """
        Makes scatter plot and saves it into png file
        """
        fig = Figure()
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_subplot('111')
        ax.set_title(figure_title)
        ax.set_xlabel(x_title)
        ax.set_ylabel(y_title)
        color_map = ['green'] * len(x_vals)
        
        xmin = min(x_vals)
        xmax = max(x_vals)
        ymin = min(y_vals)
        ymax = max(y_vals)

        if plot_std_x_axis:
            ax = self.plotThreeSTDLimits(ax, y_vals, xmin, xmax, plot_axis = "x")
        if plot_std_y_axis:
            ax = self.plotThreeSTDLimits(ax, x_vals, ymin, ymax, plot_axis = "y", plot_lower = False)
         
        if x_std_limit != None:
            mean = numpy.mean(x_vals)
            std = numpy.std(x_vals, ddof=1)
            x_limit = mean + (std * x_std_limit)
            ax.plot([x_limit, x_limit], [ymin, ymax], 'b--')
            color_map = self.getColorMap(x_vals, x_std_limit, color_map, limit_type = "larger", use_std_limit = True)
        if y_std_limit != None:
            mean = numpy.mean(y_vals)
            std = numpy.std(y_vals, ddof=1)
            y_limit = mean - (std * y_std_limit)
            ax.plot([xmin, xmax], [y_limit, y_limit], 'b--')
            color_map = self.getColorMap(y_vals, y_std_limit, color_map, use_std_limit = True)
            
        if x_val_limit != None:
            ax.plot([x_val_limit, x_val_limit], [min(y_vals), max(y_vals)], 'b--')
            color_map = self.getColorMap(x_vals, x_val_limit, color_map, limit_type = "larger")
        if y_val_limit != None:
            ax.plot([min(x_vals), max(x_vals)], [y_val_limit, y_val_limit], 'b--')
            color_map = self.getColorMap(y_vals, y_val_limit, color_map)
            
        if x_log_scale:
            ax.set_xscale('log')
        if y_log_scale:
            ax.set_yscale('log')
        ax.scatter(x_vals, y_vals, c=color_map)
        canvas.print_figure(plot_fn)
        return color_map
    
    
    
    def getColorMap(self, values, limit, col_map, limit_type = "smaller", use_std_limit = False):
        """
        Makes color map based on given limits. If standard deviation based limit is given, new limit is calculated,
        otherwise the value given is used directly. Any values smaller / larger than the limit is marked red.
        """
        if limit_type not in ("smaller", "larger"):
            raise ValueError(limit_type)
        # calculate limit if std dev is used
        if use_std_limit:
            #print ("Changing color map limit to std limit. Limit: %3.3f" % limit)
            mean = numpy.mean(values)
            std = numpy.std(values, ddof=1)
            if limit_type == "larger":
                limit = mean + (limit * std)
            else: 
                limit = mean - (limit * std)
            #print ("New std limit: %3.3f" % limit)
        #print ("color map values limit: %3.3f mean: %3.3f, min: %3.3f, max: %3.3f" % (limit, numpy.mean(values), 
        #                                                                              min(values), max(values)))
        for i, val in enumerate(values):
            if limit_type == "larger":
                if val > limit:
                    col_map[i] = "red"
            elif limit_type == "smaller":
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
        