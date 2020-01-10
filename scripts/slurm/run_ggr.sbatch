#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=ggr
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=ggr.%j.out
#################
# a file for errors from the job
#SBATCH --error=ggr.%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=02-00:00:00
#################
# Quality of Service (QOS); think of it as sending your job into a special queue; --qos=long for with a max job length of 7 days.
# uncomment ##SBATCH --qos=long if you want your job to run longer than 48 hours, which is the default for normal partition,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos, also change to normal partition
# since dev max run time is 2 hours.
##SBATCH --qos=long
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, bigmem (jobs requiring >64Gigs RAM)
#SBATCH -p akundaje,khavari
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --nodes=1
#SBATCH --exclusive
##SBATCH --ntasks-per-node=8
#################
# --mem is memory per node; default is 4000 MB per CPU, remember to ask for enough mem to match your CPU request, since
# sherlock automatically allocates 4 Gigs of RAM/CPU, if you ask for 8 CPUs you will get 32 Gigs of RAM, so either
# leave --mem commented out or request >= to the RAM needed for your CPU request.  It will also accept mem. in units, ie "--mem=4G"
##SBATCH --mem=8000
# to request multiple threads/CPUs use the -c option, on Sherlock we use 1 thread/CPU, 16 CPUs on each normal compute node 4Gigs RAM per CPU.  Here we will request just 1.
#SBATCH -c 16
#################
# Have SLURM send you an email when the job ends or fails, careful, the email could end up in your clutter folder
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=dskim89@stanford.edu
#now run normal batch commands
# note the "CMD BATCH is an R specific command

# NOTES: can only use 1 node since multiprocessing will only be on 1 node. ask for 16 even though only going to do 12 - likely gives you 20 CPUs

# cudnn and tensorflow not actually used, but potentially needed for environment
##SBATCH --gres gpu:1

# if doing anything with tronn
#module load py-tensorflow/1.7.0_py27

# plotting
module load R/3.4.0 

# for GSEA
module load java/1.8.0_131

# for IDR
module load python/3.6.1
#module load py-numpy/1.14.3_py36
#module load viz
#module load py-matplotlib/2.2.2_py36

source ~/.bashrc
source activate ggr_env

# tmp dir
TMP_DIR=$L_SCRATCH

# work dir
WORK_DIR=$PI_SCRATCH/users/dskim89/ggr/v1.1.0

# figure out a way to control total threads at top level
echo $SLURM_JOB_CPUS_PER_NODE

# actual command
ggr --cluster sherlock --threads $SLURM_JOB_CPUS_PER_NODE --out_dir $TMP_DIR
rsync -avz --progress $TMP_DIR/ $WORK_DIR/
