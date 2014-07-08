ppi.m is an implementation of the DUNE (DUplication & NEofunctional-
ization) model of eukaryotic protein-protein interaction network 
evolution.  Please see http://ppi.tinybike.net for project details.  To 
run, this script requires several supporting files:

ppi_import.m (data import, requires adjacency matrix as input file)
pointmutation.m (neofunctionalization/assimilation events)
timedepstats.m (updates dynamical properties at specified intervals)
trackstats.m (updates degree and betweenness "tracking" values)
 
In addition to these, the MatlabBGL package:
http://www.mathworks.com/matlabcentral/fileexchange/10922
and Louvain modularity algorithm:
https://sites.google.com/a/brain-connectivity-toolbox.net/bct/
are also required.

The following files can be used to calculate various statistics and plot
several static and dynamic quantities obtained from the simulation, but
are not required to run ppi.m:

std_dyn.m (calculates dynamical statistics)
std_stat.m (calculates static end-state statistics)
std_tracking.m (calculates degree and betweenness "tracking" statistics)
ppi_dyn.m (prepares dynamical quantities for plotting)
RT_save_figures_bvk.m (betweenness vs degree)
RT_save_figures_nvk.m (nearest-neighbor degree vs degree)
RT_save_figures_Cvk.m (clustering coefficients vs degree)
RT_save_figures_E2.m (second-largest eigenvalue evolution)
RT_save_figures_px.m (closeness distribution)
RT_save_figures_plambda.m (walk matrix eigenvalue distribution)
RT_save_figures_pb.m (betweenness distribution)
RT_save_figures_pk.m (degree distribution)
RT_save_figures_gcc.m (global clustering coefficient evolution)
RT_save_figures_Q.m (modularity coefficient evolution)
RT_save_figures_d.m (diameter evolution)
RT_save_figures_f1.m (largest component evolution)

Jack Peterson (jack@tinybike.net)
License: LGPL
