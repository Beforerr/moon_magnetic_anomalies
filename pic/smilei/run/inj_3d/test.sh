sim_dir=01
mkdir -p $sim_dir 
cp parameters.py test.py post_processing.ipynb _utils.py $sim_dir
cd $sim_dir
smilei test.py > job.log 2> job.err &