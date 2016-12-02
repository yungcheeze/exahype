#!/usr/bin/env python
#
# Run convergence tests.
#
#

#### COPIED FROM SHU VORTEX.
#### SHOULD BE BOTH GENERALIZED AND MERGED


"""
Run convergence tests: The EulerFlow ShuVortex.
Will setup environment variables to execute a templated specification file.
If called without arguments, it will run a series of runs in parallel.

Sample usages:

./runShuVortex.py -p 2 -m 0.1
./runShuVortex.py --all --wait
for p in 2 3 4; do ./runShuVortex -p $p -m 0.1 & done
"""

from numpy import arange
import subprocess, os, sys, time # batteries
shell = lambda cmd: subprocess.check_output(cmd, shell=True).strip()
getenv = lambda key, default=None: os.environ[key] if key in os.environ else default
me = os.path.basename(__file__)
def log(text): print "%s: %s" % (me, text)
pidlist = lambda processes: " ".join([str(proc.pid) for proc in processes])

# 3D MHD Wave
# with default MHD application but in 3D
polyorders = arange(2,10)
width = 1.0
depths = arange(1,6)

meshsize = lambda width, depth: width / 3.**depth  # actual meshsize we will get (dx)
numcells = lambda width, meshsize: width/meshsize
maxmeshsizefactor = 1.1 # to ensure meshsize is little bit smaller
maxmeshsize = lambda meshsize: meshsize * maxmeshsizefactor

settings = {}

settings['SIMBASE'] = getenv('SIMBASE', default="simulations/")
settings['ExaBinary'] = shell("echo $(exa root)/$(exa find-binary MHD)")
settings['ExaSpecfile'] = "MHD_AlfenWave3DConvergence.exahype"
settings['ExaRunner'] = "../../RunScripts/runTemplatedSpecfile.sh"

# set initial data to use.
settings['EXAHYPE_INITIALDATA']="AlfenWave"
# parameters for setting up the specfile
settings['ExaWidth']=str(width)
settings['ExaEndTime']="12.0" # MHD LONG RUN
# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
settings['ExaConvOutputRepeat']="0.1"
settings['ExaVtkOutputRepeat']="0.5"
# single threaded in the moment.
settings['ExaTbbCores'] = "1"

settings['ExaVtkFormat'] = "Legendre::vertices"

# this is useful if you compiled with assertions
settings['EXAHYPE_SKIP_TESTS'] = "True"

# template to set up a queueing system
settings['QRUNTPL'] = "srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=MHD3D-p{ExapOrder}-m{ExaMeshSize}-ConvergenceStudies"
#settings['QRUNTPL'] = ""
settings['SIMBASE'] = 'simulations/'

def start(polyorder, maxmeshsize, meshsize):
	"""
	Starts a simulation

	@param polyorder is an integer
	@param row is a row in the res pdFrame
	"""
	env = os.environ.copy()
	env.update(settings)
	env['ExaMeshSize'] = str(maxmeshsize)
	env["ExaRealMeshSize"] = str(meshsize)
	env['ExapOrder'] = str(polyorder)
	env['ExaBinary'] = '{ExaBinary}-p{ExapOrder}'.format(**env)
	env['SIMDIR']="{SIMBASE}/p{ExapOrder}-meshsize{ExaMeshSize}/".format(**env)
	if 'QRUNTPL' in env:
		env['QRUN'] = env['QRUNTPL'].format(**env)
	
	if not os.path.exists(env['ExaBinary']):
		raise IOError("Failure: '{ExaBinary}' does not exist, probably you forgot to"
			"compile for all neccessary polynomial orders before?"
			"Hint: Try `exa polycompile EulerFlow`".format(**env))

	print "Starting p={ExapOrder}, maxmeshsize={ExaMeshSize} with {ExaRunner}".format(**env)
	return subprocess.Popen([env['ExaRunner']], env=env)

def runRange(polyorders=polyorders, depths=depths):
	"""
	Having polyorders, depths being two lists, run the cartesian product
	itertools.product(polyorders, depths) of simulations.
	However, do not run too fine simulations: Make a cut off at high polynomial
	orders.
	If you call this with no arguments, it takes the global defaults which give
	a reasonable set of simulations
	"""
	import pandas as pd
	res = pd.DataFrame({'depth': depths})
	res['meshsize'] = dx = meshsize(width, depths)
	res['numcells'] = numcells(width, dx)
	res['maxmeshsize'] = maxmeshsize(dx)

	print "ExaHyPE convergence Analysis"
	print "============================"

	print "Can do convergence analysis on %fx%f sized domain with the following grid props: " % (width,width)
	print res
	print "Will do all these tests for these orders of the polynomial order:"
	print polyorders
	print "However, now reducing the runs at big resolutions"
	
	# adapt the number of runs/ maximum cells for the polynomial degree
	# adaptrunrange[<polyorder>] -> [<subset of maxmeshsizes>]
	# numcells is array([   3.,    9.,   27.,   81.,  243.])
	until = lambda maxcells: res[ res['numcells'] <= maxcells ]
	adaptrunrange = { 2: until(243), 3: until(243), 4: until(243), 5: until(243),
			6: until(81), 7: until(81),
			8: until(27), 9: until(27) }

	processes = []
	for p, rows in adaptrunrange.iteritems():
		for i, row in rows.iterrows():
			# for testing:
			#proc = subprocess.Popen(["/bin/sleep", str(p)])
			# for real:
			proc = start(p, row['maxmeshsize'], row['meshsize'])
			processes.append(proc)
	return processes

def main():
	import argparse # python included
	parser = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-p', '--polyorder', type=int, help="Polynomial order")
	parser.add_argument('-m', '--meshsize', type=float, help="Maximum mesh size (dx)")
	parser.add_argument('-w', '--wait', action='store_true', help="Wait until ExaHyPE processes have finished")
	parser.add_argument('-a', '--all', action='store_true', help='Start all sensible simulations')
	args   = parser.parse_args();
	
	if(args.polyorder or args.meshsize):
		if(not args.polyorder or not args.meshsize):
			parser.error("Please specify both the polynomial order and the meshsize or none")
		p = start(args.polyorder, maxmeshsize(args.meshsize), args.meshsize)
		processes = [p]
		log("ExaHyPE runner has been started, PID = %d" % p.pid)
	elif(args.all):
		processes = runRange()
		log("%d processes have been started with PIDS:" % len(processes))
		log(pidlist(processes))
	else:
		log("Choose either --all or -p/-m combination for run")
		parser.print_help()
		sys.exit(2)
	
	if(args.wait):
		livingprocs = processes[:]
		while True:
			# filter out finished processes
			livingprocs = [proc for proc in livingprocs if proc.poll() == None]
			log("Waiting for %d processes (%s) to finish..." % (len(livingprocs),pidlist(livingprocs)))
			if not len(livingprocs):
				break
			time.sleep(1)
			
			# Todo: Should kill processes when this script is killed inside this loop.
			# or tell processes to quit when args.wait is on and parent is killed.
			
			# do this to access values
			#proc.communicate()
		log("All processes finished. Exit codes are:")
		print [proc.returncode for proc in processes]
		# now start convergence table analysis if wanted.
	else:
		log("Done.")

if __name__ == "__main__":
	main()

# script ends here.
# Watch the further output of the simulations on your computer using 'top' or 'htop'.
# Use 'multitail simulations/*/*.log' to watch what the simulations are doing.
# When finished (or even before), use 'finish-convergence-table.py' to inspect the results.
