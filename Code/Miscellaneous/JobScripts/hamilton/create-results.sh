module load slurm



sbatch --nodes=4 --ntasks-per-node=4 euler.slurm-script

sbatch --nodes=1  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=2  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=3  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=4  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=5  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=6  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=7  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=8  --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=16 --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=32 --ntasks-per-node=4 euler.slurm-script
sbatch --nodes=64 --ntasks-per-node=4 euler.slurm-script


sbatch --nodes=1  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=2  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=3  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=4  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=5  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=6  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=7  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=8  --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=16 --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=32 --ntasks-per-node=8 euler.slurm-script
sbatch --nodes=64 --ntasks-per-node=8 euler.slurm-script

sbatch --nodes=1  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=2  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=3  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=4  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=5  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=6  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=7  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=8  --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=16 --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=32 --ntasks-per-node=16 euler.slurm-script
sbatch --nodes=64 --ntasks-per-node=16 euler.slurm-script

