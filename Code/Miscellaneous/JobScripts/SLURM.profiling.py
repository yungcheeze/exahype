from SLURM import runMultipleSlURMjobs

#processes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 1680]
#pdegrees = [1, 2, 3, 4, 5, 6, 7, 8]
#hmaxs = ["0.4", "0.2", "0.04", "0.02", "0.005", "0.002", "0.0005"]
#compilers = ["GNU", "Intel"]
#dimensions = ["2D", "3D"]
#modes = ["Profile", "Release"]

processes = [2, 64, 128]
threads = [1]
h_p_ts = [["0.04", 3, 0.001], ["0.04", 8, 0.001]]
compilers = ["GNU", "Intel"]
dimensions = ["2D"]
modes = ["Profile"]

runMultipleSlURMjobs(dimensions, processes, threads, h_p_ts, compilers, modes, "profilingRun")