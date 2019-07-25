#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_hhuang_1
# $ -cwd = run from this current working directory (relative paths)
# -o stdout-Simulation_1_rho.out
# $ -l h_data=16G,h_rt=24:00:00
# $ -t 1-10:1

source ~/.bash_profile
#SGE_TASK_ID=1
num_of_snps=(1 2 3 4 5)
this_num_of_snps=${num_of_snps[$SGE_TASK_ID-1]}

simulation='Simulation_3_numOfSnps='$this_num_of_snps
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_numOfSnps_'$SGE_TASK_ID'_results'

for i in {1..1000};do
	python3 ../simulate.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt -c $num_of_snps
	python3 ../MCaviar.py -l LD -z Z -o $simulation'/'$outfile'_'$i -c $num_of_snps
	rm -r LD
	rm -r Z

	speculate_file=$outfile'_'$i'_set.txt'
	python3 ../capture.py -s $speculate_file -t true_causal_set.txt -p '../'$simulation
done
