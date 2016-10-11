#   Request 1 processors on 1 node 
#   
#PBS -m e
#PBS -M hongyuz@andrew.cmu.edu
#
#PBS -l nodes=1:ppn=1
#
#PBS -l walltime=100:00:00
#
#   Request 1 gigabyte of memory per process
#
#   Request that regular output and terminal output go to the same file
#
#PBS -j oe
#
#   The following is the body of the script. By default,
#   PBS scripts execute in your home directory, not the
#   directory from which they were submitted. The following
#   line places you in the directory from which the job
#   was submitted.
#
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=16
export I_MPI_PIN=off

#
#
echo " "
echo " "
echo "Job started on `hostname` at `date`"
python nbody_mean.py /physics2/hongyuz/results/a_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/a_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/b_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/b_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/c_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/c_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/e_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/e_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/f_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/f_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/h_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/h_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/i_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/i_1024_5_pot_500.dat
python nbody_mean.py /physics2/hongyuz/results/j_1024_5.dat 100 100 100 100 100 100 1000 0 1 /physics2/hongyuz/results/j_1024_5_pot_500.dat

echo " "
echo "Job Ended at `date`"
echo " "
