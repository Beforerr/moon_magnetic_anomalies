export OMP_NUM_THREADS=12
export OMP_SCHEDULE=dynamic
export OMP_PROC_BIND=true
spack load miniconda3 hdf5%oneapi
mkdir "log" "result" -p
mpirun -np 1 $HOME/Smilei/smilei test.py > log/sim.log
# cp * result/1