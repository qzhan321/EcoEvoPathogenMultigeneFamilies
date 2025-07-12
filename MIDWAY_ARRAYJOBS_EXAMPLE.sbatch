#!/bin/bash
#SBATCH --job-name=Rdivtest
#SBATCH --time=36:00:00
#SBATCH --output=Rdivtest_%A_%a.out
#SBATCH --error=Rdivtest_%A_%a.err
#SBATCH --array=1-2
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30000
#SBATCH --partition=broadwl
#SBATCH --account=covid-19
#SBATCH --qos=covid-19
#SBATCH --mail-type=ALL
#SBATCH --mail-user=heqixin@uchicago.edu

# Print this sub-job's task ID
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
module load gcc/6.1
module load python/3.7.0
module load R

filePrefix=Rdivtest

#here define where the actual run will be stored
remoteDir=/scratch/midway2/heqixin/${filePrefix}_${SLURM_ARRAY_TASK_ID}

#copy original Rdiv code to the running folder
cp -r /home/heqixin/varyRecomb/varEvol $remoteDir

cd $remoteDir
python writeParameters.py -p RdivTestParam.csv -i RdivTestTemplate.py -n $SLURM_ARRAY_TASK_ID -x $filePrefix

#build the model
./build.py -p ${filePrefix}_${SLURM_ARRAY_TASK_ID}_input.py -d runfolder
#execute the run
cd runfolder
./bin/varRdiv


