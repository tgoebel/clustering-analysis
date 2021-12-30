# clustering-analysis
Seismicity Clustering Analysis Based on nearest neighbor distances of event pairs

To separate seismicity into background and clustered events, we use the distribution of nearest-neighbor event pairs and compare observed clustering characteristics with expectations from random poissonian earthquakes. Commonly, space-time distances can be described by a 2D bi-modal distribution. The first mode at small interevent times and distances highlights clustered events (e.g. aftershocks), whereas the second mode at larger distances is comprised of background events. The background mode for the California catalog corresponds to that expectation (see plots). We use the 99th percentile of nearest-neighbor distances from the randomized catalogs to separate background and clustered events, which allows for a clear separation between the two modes in California.


Dependencies:

Python 3.7
Numpy, matplotlib, matplotlib-Basemap, scipy, scipy, datetime, calendar

Use the following references if you use this code:
- Zaliapin, I., and Ben-Zion, Y., 2013, Earthquake clusters in southern California I: Identification and stability: Journal of Geophysical Research: Solid Earth, v. 118, no. 6, p. 2847–2864, doi: 10.1002/jgrb.50179.

- Goebel, T.H.W., Rosson, Z., Brodsky, E.E., and Walter, J.I., 2019, Aftershock deficiency of induced earthquake sequences during rapid mitigation efforts in Oklahoma: Earth and Planetary Science Letters, v. 522, p. 135–143, doi: 10.1016/j.epsl.2019.06.036.


# Tutorial:

1) Download standard or relocated catalog (e.g. https://service.scedc.caltech.edu/eq-catalogs/date_mag_loc.php
or https://scedc.caltech.edu/research-tools/altcatalogs.html)
2) Convert catalog to EqCat object with attribute self.data, which is a dictionary with data columns
'Time' = Decimal Year, 'Lon', 'Lat', 'Mag', 'Depth'. Use '1_create_mat_eqCat_file.py' to do the conversion and
save the EqCat as matlab binary (.mat). Alternatively, earthquake catalog formats can be changed in matlab
and saved as .mat with variable names: 'Time' = Decimal Year, 'Lon', 'Lat', 'Mag', 'Depth'.
An example catalog is provided in the /data directory (hs_1983_2011_all.mat)

The following steps require estimates of fractal dimension, D, completeness magnitude, Mc, and b-value.
Mc and the b-value can be estimated using the github repository: https://github.com/tgoebel/magnitude-distribution
It is recommended that the sensitivity of the results to changes in these parameters are tested.

3) Compute separation between clustered and background events: '2_eta_0.py'

4) Compute nearest neighbor distances (NND) and find parent event for each except for the first event in the catalog:
'3_NND.py'

5) Assemble event families and save them within cluster. Each cluster has a unique ID which is used as variable name
in a corresponding python dictionary. The clusters contain the unique event IDs of all members. Note that the cluster
with ID and variable name '0' contains singles, i.e. events with parents at nearest-neighbor-distance beyond eta_0:
'6_createClust.py'

6) Count the numbe of events within each cluster (or family) before (foreshocks) and after (aftershocks) the largest
magnitude event in each family. Singles have 0 fore and aftershocks:
'7_productivity.py'

7) Plot productivity relationship including number of aftershocks in ech family and average number of aftershocks
within magnitude-binned mainshocks. Plot alpha=1 slope for comparison:
'8_plot_productivity.py'

All results should be compared to the provided figures using the scripts:
'1b_plot_eqLocs.py', '4_dist_tau.py', '5_plot_lat_t.py'
- also check results in:
  Zaliapin, I., and Ben-Zion, Y., 2013, Earthquake clusters in southern California I: 
  Identification and stability: Journal of Geophysical Research: Solid Earth, v. 118, no. 6, p. 2847–2864, doi: 10.1002/jgrb.50179.
  and 
- Goebel, T.H.W., Rosson, Z., Brodsky, E.E., and Walter, J.I., 2019, Aftershock deficiency of induced earthquake sequences during rapid mitigation efforts in Oklahoma: Earth and Planetary Science Letters, v. 522, p. 135–143, doi: 10.1016/j.epsl.2019.06.036.
