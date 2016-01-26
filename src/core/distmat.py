'''
Created on 26.1.2016

@author: paakkone
'''
import os
import numpy

class DistanceMatrix(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        pass
    
    def calcDistMat(self, x, y):
        print ("Max x-axis: %s, max y-axis: %s" % (max(x), max(y)))
        x = numpy.array(x)
        y = numpy.array(y)
        x_range = max(x) - min(x)
        y_range = max(y) - min(y)
        x = (x - min(x)) / x_range
        y = (y - min(y)) / y_range
        print ("After normalization: Max x-axis: %s, min %s,  max y-axis: %s, min %s" % (max(x), min(x), max(y), min(x)))
        dx = x[..., numpy.newaxis] - x[numpy.newaxis, ...]
        dy = y[..., numpy.newaxis] - y[numpy.newaxis, ...]
        d = numpy.array([dx, dy])
        return (d**2).sum(axis=0)**0.5
    
    def calcDistMat2(self, x, y):
        vals = []
        # normalize
        max_x = max(x)
        max_y = max(y)
#         for i in range(len(x)):
#             vals.append((x[i],y[i]))
        for i in range(len(x)):
            vals.append((x[i]/max_x,y[i]/max_y))
        matrix = numpy.array([[complex(val[0], val[1]) for val in vals]])
        dist_mat = abs(matrix.T - matrix)
        return dist_mat
    
    def getDistanceMatrix(self, x, y):
        """
        Make distance matrix and get median distance
        """
        print ("len(x): %d, len(y): %d" % (len(x), len(y)))
        dist_mat = self.calcDistMat(x,y)
        mean_val = dist_mat.mean()
        std_val = numpy.std(dist_mat)
        median_val = numpy.median(dist_mat)
        print ("Minimum distance: %f, maximum_distance: %f, mean: %f, median: %f" % (dist_mat.min(), dist_mat.max(), mean_val, median_val))
        print ("STD: %s" % std_val)
        no_short_distances = numpy.zeros(len(dist_mat[0]))
        #mad_x = numpy.zeros(len(x))
        #mad_y = numpy.zeros(len(y))
        total_sums = numpy.zeros(len(dist_mat[0]))
        #print ("Dist matrix:\n%s" % dist_mat)
        #print ("dist matrix [0]\n%s" % dist_mat[0])
        #print ("len first row: %d" % len(dist_mat[0]))
        for i in range(len(x)):
            #no_short_distances[i] = numpy.sum(dist_mat[i] < (mean_val + (0 * std_val)))
            no_short_distances[i] = numpy.sum(dist_mat[i] < (mean_val / 1 ))
            total_sums[i] = numpy.sum(dist_mat[i])
            #print ("Row: %d, COUNT: %s" % (i, count))
            #print ("len(count): %d, len(count[0]): %d, sum: %s" % (len(count), len(count[0]), sum(count)))
        ave_short = numpy.average(no_short_distances)
        std_short = numpy.std(no_short_distances)
        std_mult = 3.0
        limit = ave_short - (std_mult * std_short)
        count2 = numpy.sum(no_short_distances < limit)
        print ("AVE: %f, STD: %f, Count2: %s" % (ave_short, std_short, count2))
        points = no_short_distances < limit
        print ("len(points): %d, points: %s, sum: %s, limit: %s" % (len(points), points, numpy.sum(points), limit))
        
#         total_sum_ave = numpy.average(total_sums)
#         total_sum_std = numpy.std(total_sums)
#         print("Total sum min: %s, max: %s" % (numpy.min(total_sums),numpy.max(total_sums)))
#         print ("Total sum: AVE: %f, STD: %f" % (total_sum_ave, total_sum_std))
#         limit = total_sum_ave + (std_mult * total_sum_std)
#         points = total_sums > limit
#         print ("len(points): %d, points: %s, sum:%s, limit:%s" % (len(points), points, numpy.sum(points), limit))

        colors = []
        col1= "green"
        col2 = "red"
        for p in points:
            if p == True:
                colors.append(col2)
            else:
                colors.append(col1)
        return colors
