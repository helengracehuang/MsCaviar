#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_hhuang_1
# $ -cwd = run from this current working directory (relative paths)
# $ -o stdout-Simulation_1_numOfStudies.out
# $ -l h_data=4G,h_rt=24:00:00
# $ -t 1-5:1

source ~/.bash_profile
# SGE_TASK_ID=1
num_of_studies=(1 2 3 4 5)
this_num=${num_of_studies[$SGE_TASK_ID-1]}

simulation='Simulation_2_numOfStudies='$this_num
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_numOfStudies_'$SGE_TASK_ID'_results'

for i in {1..1000};do
	python3 ../simulate.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt -n $this_num
	python3 ../MCaviar.py -l LD -z Z -o $outfile'_'$i
	rm -r LD
	rm -r Z

	speculate_file=$outfile'_'$i'_set.txt'
	python3 ../capture.py -s $speculate_file -t true_causal_set.txt -p '../'$simulation
done
