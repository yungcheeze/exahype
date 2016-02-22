1. Installation:
Copy the bash scripts in the application folder

2. Running the scripts:
On the command line, execute:

source config_<project>
BASH_ENV=params_<system_id> ./job_scheduler.sh

Example for euler2d on Haswell: 
source config_eulerflow2d
source params_haswell
./job_scheduler.sh

The job scheduler will create a subfolder
scaling in the project folder (eulerflow2d in the example above),
and further a subfolder with a name consisting
of the current date and the system id.
Within the last subfolder, the script stores
the ExaHyPE output files.
