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
        
    def plotSTDLimits(self, ax, x_values, y_values, x_val_limit, no_std, line_color, color_map = True):
        mean_heterozygosity = numpy.mean(y_values)
        std_heterozygosity = numpy.std(y_values, ddof=1)
        upper = mean_heterozygosity + (no_std * std_heterozygosity)
        lower = mean_heterozygosity - (no_std * std_heterozygosity)
        ax.plot([min(x_values), max(x_values)], [upper,upper], line_color)
        ax.plot([min(x_values), max(x_values)], [lower,lower], line_color)
        col_map = []
        if color_map:
            for i in range(len(y_values)):
                y_val = y_values[i]
                x_val = x_values[i]
                if x_val_limit == None:
                    if y_val > upper or y_val < lower:
                        col_map.append("red")
                    else:
                        col_map.append("green")
                else:
                    if y_val > upper or y_val < lower or x_val > x_val_limit:
                        col_map.append("red")
                    else:
                        col_map.append("green")
        return ax, col_map
        
    def plotScatterPlot(self, plot_fn, x_vals, y_vals, x_val_limit = None, y_val_limit = None, chosen_std_limit = 4, color_map = None, figure_title = "Title", 
                        x_title = "x-axis", y_title = "y_axis", x_log_scale = True, y_log_scale = False):
        """
        Makes scatter plot and saves it into png file
        """
        fig = Figure()
        canvas = FigureCanvasAgg(fig)
        ax = fig.add_subplot('111')
        ax.set_title(figure_title)
        ax.set_xlabel(x_title)
        ax.set_ylabel(y_title)
        
        color_map = None
        stds_and_colors = ((5, "g--"), (4, "y--"), (3, "r--"))
        for std, col in stds_and_colors:
            ax, col_map = self.plotSTDLimits(ax, x_vals, y_vals, x_val_limit, std, col)
            if std == chosen_std_limit:
                color_map = col_map
                
        if x_val_limit != None:
            ax.plot([x_val_limit,x_val_limit], [min(y_vals),max(y_vals)], 'b--')
            if color_map == None:
                color_map = self.getColorMap(x_vals, x_val_limit)
        if y_val_limit != None:
            ax.plot([min(x_vals),max(x_vals)], [y_val_limit,y_val_limit], 'b--')
            if color_map == None:
                color_map = self.getColorMap(y_vals, y_val_limit)
        if x_log_scale:
            ax.set_xscale('log')
        if y_log_scale:
            ax.set_yscale('log')
        ax.scatter(x_vals, y_vals, c=color_map)
        canvas.print_figure(plot_fn)
        return color_map
    
    def getColorMap(self, values, limit, limit_type = "smaller"):
        col_map = []
        for val in values:
            if limit_type == "larger":
                if val > limit:
                    col_map.append("red")
                else:
                    col_map.append("green")
            elif limit_type == "smaller":
                if val < limit:
                    col_map.append("red")
                else:
                    col_map.append("green")
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
        