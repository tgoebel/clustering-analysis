# python2.7
"""
        - functions to plot binned and smoothed desnities of 2D data points
        
            --> plot as 2D probability density with Sum^x_y = 1, Integral = 1
            -For that purpose: divide by current sum (~number of events) and dx, dy
               which are the bins in x and y
"""


import numpy as np
import matplotlib.pyplot as plt

np.random.seed(12345)
#================================================================================
#                           fct. definitions
#================================================================================   
def density_2D( x, y, x_bin, y_bin, **kwargs):
    """
        2D, smoothed event density for point cloud with coordinates x,y
        uses method: scipy.stats.kde.gaussian_kde
    :input     x,y            - dataset
               x_bin, y_bin   - binned x and y vectors
    
    
        kwargs['sigma'] - specify gaussian smoothing kernel ('bw_method' in scipy.stats.kde)
                          default: = 'scott' 
                              sigma = n**( -1./(d+4)), d- number of dimensions, n - number of data points
                        - 'silverman'
                              sigma = (n * (d + 2) / 4.)**(-1. / (d + 4))
                        - float( )   = set Gaussian Bandwidth directlty
                                        
                           
    return XX, YY, ZZ - 2D binned x and y coordinates and density for each cell
    """
    from scipy.stats import kde
    sigma = 'scott'
    if 'sigma' in kwargs.keys() and kwargs['sigma'] is not None:
        sigma = kwargs['sigma']
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    fct_Gauss2D = kde.gaussian_kde( np.array([x,y]), bw_method = sigma)
    # meshgrid of x and y coordinates
    XX,YY   = np.meshgrid( x_bin, y_bin)
    ZZ = fct_Gauss2D( np.vstack([XX.flatten(), YY.flatten()])).reshape( XX.shape)
    dx, dy = x_bin[1] - x_bin[0], y_bin[1] - y_bin[0]
    # check if integral is ~ zero, better: use midepoint method
    print( 'check if integral ~1', ZZ.sum()*( dx*dy)) #ZZ[ZZ>0].mean()*(XX.max()-XX.min())*(YY.max()-YY.min()))
    return XX-.5*dx, YY-.5*dy, ZZ
    #return XX, YY, ZZ
#================================================================================
#                          parameters
#================================================================================  
Nev   = 10
nbins = 30
xmin, xmax = -2.5, 2.5
ymin, ymax = -5,    5

sigma = .1 # Gaussian smoothing kernel

# binsize in x and y
dx, dy= float(xmax-xmin)/nbins, float(ymax-ymin)/nbins

#================================================================================
#                    create random data, and binned data
#================================================================================   
# Create data: 200 points
data = np.random.multivariate_normal([0, 0], [[1, 0.5], [0.5, 3]], Nev)
x, y = data.T
a_xbin = np.arange( xmin-dx, xmax+2*dx, dx) 
a_ybin = np.arange( ymin-dy, ymax+2*dy, dy)

#================================================================================
#                     compute Gaussian density
#================================================================================ 
XX,YY,ZZ = density_2D( x, y, a_xbin, a_ybin, sigma = sigma)
#================================================================================
#                           plots
#================================================================================ 
# Create a figure with 6 plot areas
fig, axes = plt.subplots(ncols=6, nrows=1, figsize=(21, 5))
 
# Everything sarts with a Scatterplot
axes[0].set_title('Scatterplot')
axes[0].plot(x, y, 'ko')
# As you can see there is a lot of overplottin here!
 
# Thus we can cut the plotting window in several hexbins

axes[1].set_title('Hexbin')
axes[1].hexbin(x, y, gridsize=nbins, cmap=plt.cm.BuGn_r)
axes[1].plot( x, y, 'ko', ms= 2)
axes[1].set_xlim( axes[0].get_xlim())
axes[1].set_ylim( axes[0].get_ylim())

# 2D Histogram
axes[2].set_title('2D Histogram')
counts, xedges, yedges, __ = axes[2].hist2d( x, y, bins=nbins, cmap=plt.cm.BuGn_r, normed = True)
axes[2].plot( x, y, 'ko', ms= 2)
axes[2].set_xlim( axes[0].get_xlim())
axes[2].set_ylim( axes[0].get_ylim())

dx,dy = (xedges[1]-xedges[0]), (yedges[1]-yedges[0])
#print(xedges, yedges, counts)
print( 'check if integral ~1', counts.sum()*( dx*dy), counts.mean()*(xedges[-1]-xedges[0])*(yedges[-1]-yedges[0]))


 
# plot a density
axes[3].set_title('Gaussian Smoothing')
axes[3].pcolormesh( XX, YY, ZZ, cmap=plt.cm.BuGn_r)
axes[3].plot( x, y, 'ko', ms= 2)
axes[3].set_xlim( axes[0].get_xlim())
axes[3].set_ylim( axes[0].get_ylim())



# add shading
axes[4].set_title('2D Density with shading')
axes[4].pcolormesh( XX,YY,ZZ, shading='gouraud', cmap=plt.cm.BuGn_r)
axes[4].set_xlim( axes[0].get_xlim())
axes[4].set_ylim( axes[0].get_ylim())

 
# contour
axes[5].set_title('Contour')
axes[5].pcolormesh( XX, YY, ZZ, shading='gouraud', cmap=plt.cm.BuGn_r)
axes[5].contour(     XX, YY, ZZ )
axes[5].set_xlim( axes[0].get_xlim())
axes[5].set_ylim( axes[0].get_ylim())
plt.show()










