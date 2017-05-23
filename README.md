# tracking-and-frap
A set of tools for analyzing diffusion and FRAP data.

After discovering the incredible TrackMate plugin for ImageJ (https://github.com/fiji/TrackMate)
I've decided to write a few modules for parsing its outputs in Python and analyzing/plotting them 
in ways that are not included in the plugin itself. The list will hopefully keep growing, but 
here's what I want this project to be able to do for now:

1) Parse TrackMate's output xml files and reorganize the data into track-based lists of spot
coordinates (provided the tracks don't branch). This is useful, e.g., for diffusion analysis.
2) Allow for custom MSD analysis, e.g. grouping particles by intensity or radius to see how
their movement changes in response to these properties.
3) Simplify analysis of FRAP experiments of moving spots.
4) and hopefully much more, TBD
