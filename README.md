# Cluster analysis
A set of tools for analyzing the distribution and dynamics of molecular clusters.

Many components of living cells, such as signaling proteins, form dynamic clusters or oligomers
in response to changes in the cellular environment. These clusters can typically be visualized
as bright puncta by fluorescence microscopy. Unfortunately, more often than not such clusters 
are reported and analyzed in a superficial and descriptive way that misses the exciting 
underlying biophysical processes that lead to their formation, trafficking, and dissolution.

This project was started to create a set of convenient tools and analysis approaches that
can be used to objectively and reproducibly quantify a wide range of spatial and dynamic parameters
of molecular clusters. For now, I'm mainly just using this for my research and am adding new code
whenever I need to answer a new specific question about a real-world dataset, but since many of the
approaches are readily generalizable to any clustering molecules, I believe some of the code may be 
useful to other scientists in the long run.

Most of the analysis code is written in Python, but the project also includes short ImageJ scripts
and CellProfiler pipelines to handle the first few steps of image processing.

Below is the ever-growing list of functionality available in the project:

1) Parsing the output xml files of the excellent ImageJ plugin TrackMate 
(https://github.com/fiji/TrackMate) output xml files and reorganize the data into track-based 
lists of spot coordinates (provided the tracks don't branch). This is useful, e.g., for diffusion analysis.

2) Custom MSD analysis on particle trajectories, e.g. grouping particles by intensity or radius to see how
their movement changes in response to these properties.

3) Automated FRAP analysis on small diffusing spots with local background correction.

4) General FRAP analysis.

5) Plotting the distribution of clusters per cell in high-content imaging experiments. This
makes heavy use of the free and open source CellProfiler suite (https://cellprofiler.org/)
for the ground work of identifying cells and extracting cellular features.

6) A set of tools for parsing the raw HDF5 output from CellProfiler runs.

7) And hopefully much more to come!
