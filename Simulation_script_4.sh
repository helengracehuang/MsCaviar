#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_rosemary_1
# $ -cwd = run from this current working directory (relative paths)
# -o stdout-Simulation_4_tau_sqr.out
# $ -l h_data=6G,h_rt=24:00:00
# $ -t 1-9:1

source ~/.bash_profile
#SGE_TASK_ID=1
tau_2=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
this_tau_2=${tau_2[$SGE_TASK_ID-1]}

simulation='Simulation_4_tau_sqr='$this_tau_2
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_tau_sqr_'$SGE_TASK_ID'_results'

for i in {1..1000};do
	python3 ../simulate.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt -t $this_tau_2
	python3 ../MCaviar.py -l LD -z Z -o $outfile'_'$i -t $this_tau_2
	rm -r LD
	rm -r Z

	speculate_file=$outfile'_'$i'_set.txt'
	python3 ../capture.py -s $speculate_file -t true_causal_set.txt -p '../'$simulation
done
