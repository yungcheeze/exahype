# The Exa Python Postprocessing Library

This directory holds various Python files helpful for postprocessing of
ExaHyPE VTK simulation results. These tools allow stripping down the
complicated VTK structure to simple CSV/ASCII files or NumPy tables,
suitable for further processing with gnuplot or Python.

This is especially interesting for 1D and 2D slices to understand what
is going on in a simulation.

The Python scripts can be used both as command line programs as well as
python libraries for interactive scripting and plotting. The name *Exa*
is not really a joke, but of course this is not high-performance.

## Usage/Dependencies

You just need a scientific python installation, ie. Anaconda or the
usual Ubuntu packages `python-numpy python-scipy python-matplotlib`.
You also need VTK for Python, ie. the package `python-vtk` at
Ubuntu.
