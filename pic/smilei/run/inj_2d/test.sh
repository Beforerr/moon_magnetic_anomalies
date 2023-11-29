sim_dir=05
mkdir -p $sim_dir 
cp parameters.py test.py post_processing.ipynb _utils.py $sim_dir
cd $sim_dir

spack env activate smilei-oneapi
export OMP_NUM_THREADS=8
export OMP_SCHEDULE=dynamic

$HOME/Smilei/build/smilei-oneapi/smilei test.py > job.log 2> job.err &