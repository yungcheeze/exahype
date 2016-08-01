<h2>File overview</h2>
The rst files are Python documentation files. The makefile is required to generate this documentation.




<h2>Available tools</h2>
Run .py -h for more details.


plot-multicore-speedup.py   Run plot-multicore-speedup.py -h for more details.
hpclib.py                   Helper routines 
conf.py                     A helper file for the documentation



<h2>Build the documentation</h2>
The Makefile is only used for building
the documentation.

To do so, run:
make docdir
make html

The first command does only
create the folder 'doc' and
three subfolders.

To delete the documentation
make clean
