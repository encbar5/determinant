#!/bin/bash
#
##################################################################
#
# GENERAL INFORMATION ABOUT THE SLURM BATCH SYSTEM ON SCHOONER
#
# greetings.sbatch
#
# last modified by Henry Neeman 2017-01-17
#
# PLEASE DON'T REMOVE THESE COMMENTS EVER!!!
# They will always be important and may save you lots of grief.
#
# Also, please note that lines that begin with pound-SBATCH (#SBATCH)
# are batch scheduler directives, so they are absolutely crucial.
# DON'T REMOVE THE pound sign (#) from before the SBATCH!!!!
#
# Everywhere throughout this batch script file:
# * CHANGE oucs001 to your username;
# * CHANGE en@ou.edu to your e-mail address.
#
# NOTE: If you create any file of human-readable text on a Windows PC,
# you *MUST* perform the following command on it:
#
#   dos2unix filename
#
# where you replace "filename" with the name of the file; for example:
#
#   dos2unix /home/oucs001/SIPE/Greetings/C/greetings_input.txt
#
# This is because virtually all text editors in Windows embed hidden
# special characters in text files (for example, font information),
# and these hidden special characters cause Unix/Linux programs to
# choke.
#
# To submit a batch job:
#
#   sbatch greetings.sbatch
#
# Note that there is no "<" sign in the above command.
#
# To see what batch jobs you personally have in the scheduler:
#
#   squeue -u oucs001
#
# To see what batch jobs everyone has in your partitions:
#
#   squeue
#
# To kill a batch job (replace JOBID with the batch job ID number):
#
#   scancel JOBID
#
# but again replacing JOBID with the batch job identifer.
#
# Below is the output of squeue.
#
# JOBID PARTITION     NAME     USER ST TIME NODES NODELIST(REASON)
#505994    normal HERV_Gag   nsloat PD 0:00    16 (Resources)
#506260    normal  wrf_No3   jdduda PD 0:00    10 (Resources)
#507599    normal may29-qn dphoenix PD 0:00    22 (Resources)
#507621    normal wrf-cbmz dphoenix PD 0:00    30 (Resources)
#511078    normal glb6x4-o scrowell PD 0:00     1 (Resources)
#511079    normal glb6x4-o scrowell PD 0:00     1 (Resources)
#511080    normal glb6x4-o scrowell PD 0:00     1 (Resources)
#512000    normal     RUN5   pattar PD 0:00     1 (Resources)
#512031    normal    RUN10   pattar PD 0:00     1 (Resources)
#513543    normal 48_prime   soumya PD 0:00     2 (Resources)
#513743    normal     pigz  lreames PD 0:00     1 (Resources)
#514814    normal  Fa2Tre1 jtcooper PD 0:00     1 (Resources)
#515761    normal     ctrl jlabriol PD 0:00    24 (Resources)
#515762    normal    shape jlabriol PD 0:00    24 (Resources)
#515764    normal  shapeXL jlabriol PD 0:00    24 (Resources)
#516026    normal maxvg_di   yinmeg PD 0:00     1 (Resources)
#529957    normal       f2 gongfeng PD 0:00     2 (Resources)
#529962    normal       f3 gongfeng PD 0:00     2 (Resources)
#531099       ieg sbatch_s   tianrm PD 0:00     1 (ReqNodeNotAvail...)
#531108       ieg sbatch_s   tianrm PD 0:00    10 (ReqNodeNotAvail...)
#531123       ieg sbatch_s   tianrm PD 0:00     1 (ReqNodeNotAvail...)
#...
#
# Note that the ST column has the status of each batch job.
#
# Specifically, PD means PENDING (sitting in the batch queue
# waiting for that batch job's turn to begin running), R means RUN
# and CG means COMPLETING (the executable has finished running, so
# the SLURM batch system is cleaning up and preparing to shut down
# the batch job.
#
##################################################################
#
# SLURM BATCH DIRECTIVES
#
# The SBATCH directive below says the name of the partition (batch
# queue) to be used. In most cases, you should use the queue named
# normal. (For some instances of this exercise, we've created a
# special queue named sipe, but that's not currently active.)
#
#SBATCH --partition=normal # sipe
#
# The SBATCH directive below says to allow only this job to run on
# this job's compute node(s), so no other jobs, even by this user,
# can run on those nodes at the same time. DON'T CHANGE THIS.
#
#SBATCH --exclusive
#
# The SBATCH directive below says to use 40 CPU cores of one CPU chip
# on one compute node, meaning that this batch jobs is non-parallel
# (serial).
#
#SBATCH --ntasks=1
#
# The SBATCH directive below says to use 20 CPU cores (on two CPU
# chips) on each of the job's compute nodes, meaning that this batch
# jobs is parallel.
#
#SBATCH --ntasks-per-node=1
#
# The SBATCH directive below says, the working directory for
# this batch job will be the directory listed after the
# equals sign.
#
# This accomplishes the equivalent of doing a cd ("change
# directory") command into that directory, except that, by
# using the directive below, then the SLURM batch system
# will know the working directory, which can be helpful in
# properly managing your batch job.
#
#SBATCH --workdir=/home/oucs001/Competition/fortran
#
# Standard output (stdout) is the output that normally would go to
# the terminal screen.
#
# Standard error (stderr) is like stdout), except that it contains
# error messages instead of regular output.
#
# The SBATCH directive below says, send stdout and stderr to the
# filenames listed below.
#
# Note that, in these filenames, %J will be replaced by the batch
# job ID number.
#
# Everywhere throughout this batch script file:
# * CHANGE oucs001 to your username;
# * CHANGE en@ou.edu to your e-mail address.
#
#SBATCH --output=/home/oucs001/Competition/fortran/output/Chol/det_%J_stdout.txt
#SBATCH --error=/home/oucs001/Competition/fortran/output/Chol/det_%J_stderr.txt
#
# The SBATCH directive below says to run for up to 30 minutes
# (that is, zero hours plus 30 minutes plus zero seconds) of
# "wall clock" time (time experienced in real life).
#
# Currently, the maximum allowable wall clock time per batch job
# is 48 hours (48:00:00).
#
# Acceptable time formats include:
#   minutes
#     example: 30 means 30 minutes
#   minutes:seconds
#     example: 24:30 means 24 minutes plus 30 seconds
#   hours:minutes:seconds
#     example: 47:30:22 means 47 hours plus 30 minutes plus 22 seconds
#   days-hours
#     example: 1-12 means 1 day plus 12 hours
#   days-hours:minutes
#     example: 1-12:30 means 1 day plus 12 hours plus 30 minutes
#   days-hours:minutes:seconds
#     example: 1-12:30:18 means 1 day plus 12 hours plus 30 minutes
#       plus 18 seconds
#
#SBATCH --time=00:10:00
#
# The SBATCH directive below says the name of the batch job, as it
# will appear in the batch partition listing when you do an squeue
# command.
#
# Everywhere throughout this batch script file:
# * CHANGE oucs001 to your username;
# * CHANGE en@ou.edu to your e-mail address.
#
#SBATCH --job-name=det
#
# The SBATCH directive below says the e-mail address to send
# notifications to, which should be changed to your e-mail address.
#
# Everywhere throughout this batch script file:
# * CHANGE oucs001 to your username;
# * CHANGE en@ou.edu to your e-mail address.
#
#SBATCH --mail-user=en@ou.edu
#
# The SBATCH directive below says to e-mail a notification when the
# batch job either completes or fails.
#
# If you only want e-mails when when a batch job fails, then in the
# SBATCH directive below, change ALL to FAIL.
#
#SBATCH --mail-type=FAIL
#
##################################################################
#
# UNIX SHELL COMMANDS

# Output the working directory, as a sanity check.

echo "Working directory:"
pwd
echo ""

# Set the runtime environment to match the compile environment.

module load OpenMPI
module load ScaLAPACK

# Run the executable, redirecting input from the given file.
# The date commands and the time command help track runtime
# and performance.
#
# NOTE: The line with a "<" sign, immediately before the
# second date command, redirects standard input
# (typically a user typing at the keyboard) from a file of
# human readable text.
#
# Everywhere throughout this batch script file:
# * CHANGE oucs001 to your username;
# * CHANGE en@ou.edu to your e-mail address.

N_PROC=1

for x in 100 500 1000 2000 3000 4000 5000; do
    mpirun -n $N_PROC --bind-to core /home/oucs001/Competition/fortran/detChol.out $x
done
