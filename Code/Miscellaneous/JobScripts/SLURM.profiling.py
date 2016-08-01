from SLURM import runMultipleSlURMjobs

#processes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 1680]
#pdegrees = [1, 2, 3, 4, 5, 6, 7, 8]
#hmaxs = ["0.4", "0.2", "0.04", "0.02", "0.005", "0.002", "0.0005"]
#compilers = ["GNU", "Intel"]
#dimensions = ["2D", "3D"]
#modes = ["Profile", "Release"]

processes = [2, 64, 128]
pdegrees = [3, 8]
hmaxs = ["0.04"]
compilers = ["GNU", "Intel"]
dimensions = ["2D", "3D"]
modes = ["Profile"]

runMultipleSlURMjobs(dimensions, processes, pdegrees, hmaxs, compilers, modes)