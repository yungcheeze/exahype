from SLURM import runMultipleSlURMjobs

processes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 1680]
threads = [1]
pdegrees = [1, 2, 3, 4, 5, 6, 7, 8]
hmaxs = ["0.04", "0.02", "0.005"]
compilers = ["GNU", "Intel"]
dimensions = ["2D", "3D"]
modes = ["Release"]


processes = [16, 32, 64, 128, 256, 512]
threads = [1]
h_p_ts = [["0.01", 3, 0.01], ["0.01", 9, 0.001] ]
compilers = ["Intel"]
dimensions = ["2D"]
modes = ["Release"]

processes = [16, 32, 64, 128, 256, 512]
threads = [1]
h_p_ts = [["0.01", 3, 0.01], ["0.01", 5, 0.001], ["0.01", 7, 0.001], ["0.01", 9, 0.001] ]
compilers = ["Intel"]
dimensions = ["2D"]
modes = ["Release"]

runMultipleSlURMjobs(dimensions, processes, threads, h_p_ts, compilers, modes, "scalingRun_with_tee")
