# clustering-analysis
Seismicity Clustering Analysis Based on nearest neighbor distances of event pairs

Dependencies:

Python2.7 (should run with python 3 after slight modifications)
Numpy, matplotlib, matplotlib-basemap, scipy, scipy, datetime, calendar

Use the following references if you use this code:
- Zaliapin, I., and Ben-Zion, Y., 2013, Earthquake clusters in southern California I: Identification and stability: Journal of Geophysical Research: Solid Earth, v. 118, no. 6, p. 2847–2864, doi: 10.1002/jgrb.50179.

- Goebel, T.H.W., Rosson, Z., Brodsky, E.E., and Walter, J.I., 2019, Aftershock deficiency of induced earthquake sequences during rapid mitigation efforts in Oklahoma: Earth and Planetary Science Letters, v. 522, p. 135–143, doi: 10.1016/j.epsl.2019.06.036.


Tutorial:

- run the provided example scripts in sequential order and compare results to provided figures (see /plots/)

- also check results in:
  Zaliapin, I., and Ben-Zion, Y., 2013, Earthquake clusters in southern California I: 
  Identification and stability: Journal of Geophysical Research: Solid Earth, v. 118, no. 6, p. 2847–2864, doi: 10.1002/jgrb.50179.
 
- convert seismicity catalog to EqCat() python object and save as .mat
- repeat analysis steps for new catalog
