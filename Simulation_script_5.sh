#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_rosemary_1
# $ -cwd = run from this current working directory (relative paths)
# -o stdout-Simulation_5_NCP.out
# $ -l h_data=6G,h_rt=24:00:00
# $ -t 1-9:1

source ~/.bash_profile
#SGE_TASK_ID=1
NCP=(3 5 7 9 11 13 15 17 19)
this_NCP=${NCP[$SGE_TASK_ID-1]}

simulation='Simulation_5_NCP='$this_NCP
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_NCP_'$SGE_TASK_ID'_results'

for i in {1..1000};do
	python3 ../simulate.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt -m $this_NCP
	python3 ../MCaviar.py -l LD -z Z -o $outfile'_'$i
	rm -r LD
	rm -r Z

	speculate_file=$outfile'_'$i'_set.txt'
	python3 ../capture.py -s $speculate_file -t true_causal_set.txt -p '../'$simulation
done
