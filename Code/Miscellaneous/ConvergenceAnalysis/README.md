Convergence Analysis for ExaHyPE
================================

These scripts organize convergence order tests in a way like they
are done by Michael Dumbser, ie. for a given polynomial order it
runs stuff for different grid sizes.

With the Shu Vortex, the table should look eg. like the following:

```
p2p2
						
numPts	l1      	l2      	linf    	o1  	o2  	oinf
20	2,43E-02	6,28E-03	5,29E-03
40	5,03E-03	1,27E-03	1,16E-03	2,27	2,30	2,19
60	1,98E-03	4,91E-04	4,41E-04	2,30	2,35	2,39
80	1,01E-03	2,46E-04	2,21E-04	2,35	2,40	2,40
```

The quantity `oI` with `I` in `{1, 2, inf}` is computed with the equation

```
  oI[row i] = log(lI[row i] / lI[row i-1]) / log(numPts[row i-1] / numPts[row i])
```

Number of points in ExaHyPE
===========================

In ExaHyPE, the number of points can only be one of (1./3.)**n, n element N,
thus possible number of points are

```
   3, 9, 27, 81, 243, ...
```

corresponding to the fixed mesh sizes

```
   0.333, 0.111,  0.03703704,  0.01234568  0.00411523, ...
```

This means that while in Michaels example, `log(numPts[i-1]/numPts[i])=log(2)`,
here this is always `log(3)`, which should not make any difference, thought.


Usage of the scripts
====================

Each script can be used for itself, but this is the overall idea:

* Start with `start-convergence.py` by adopting the parameters to your needs.
  When executing, it will start all exahype instances at the same time.
* Adapt `run-convergence.sh` to your needs, it is just a wrapper for the
  ExaHyPE binary and configuration handling.
* Use `finish-convergence-table.py` to collect the results and
  compute the convergence table in a similar fashion as done by Michael Dumbser
  in Excel.


Example of starting the code
===========================


Output of start*.py script:

```
ExaHyPE convergence Analysis
============================
Can do convergence analysis on 15.000000x15.000000 sized domain with the following grid props:
   depth  meshsize  numcells  maxmeshsize
0      1  5.000000         3     5.500000
1      2  1.666667         9     1.833333
2      3  0.555556        27     0.611111
3      4  0.185185        81     0.203704
4      5  0.061728       243     0.067901
Will do all these tests for these orders of the polynomial order:
[2 3 4 5 6 7 8 9]
Starting p=2, maxmeshsize=5.500000 
Starting p=2, maxmeshsize=1.833333 
WIPING existing simulation at simulations//p2-meshsize5.5/
Starting p=2, maxmeshsize=0.611111 
WIPING existing simulation at simulations//p2-meshsize1.83333333333/
...
```
