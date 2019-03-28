# clustering-analysis
Seismicity Clustering Analysis Based on nearest neighbor distances of event pairs

@author tgoebel - UC - Santa Cruz

This repository is an implementation of Zaliapin & Ben-Zion in python.

When using this code, please cite:
- Zaliapin, I., and Ben-Zion, Y., 2013, Earthquake clusters in southern California I: Identification and stability: Journal of Geophysical Research: Solid Earth, v. 118, no. 6, p. 2847–2864, doi: 10.1002/jgrb.50179.
and
- Goebel, T.H.W., Rosson, Z., Brodsky, E.E., and Walter, J.I., 2019, Aftershock deficiency of induced earthquake sequences during rapid mitigation efforts in Oklahoma: EPSL, v. submitted.

The provided codes reproduce Figures 4, 5 and 8 in the paper, as well as general
productivity and aftershock rate plots for Southern California.

The analysis is based on a relocated earthquake catalog from Hauksson et al. 2012 from 1981 to 2011,
available at: http://scedc.caltech.edu/research-tools/altcatalogs.html

The two files provides are: hs_1981_2011_all.mat
in matlab binary format which can be loaded using scipy.io.loadmat

and the catalog in original format:
hs_1981_2011_all.txt 
with the data columns:
#YR   MO DY HR MN SC          N     Lat         Lon      Depth   MAG  NP Dist  rms   d_n meth clID  nCl     ndTT  err_h   err_z   err_h2  err_z2 type  meth2  poly                  


