#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_hhuang_1
# $ -cwd = run from this current working directory (relative paths)
# $ -o stdout-Simulation_1_rho.out
# $ -l h_data=16G,h_rt=24:00:00
# $ -t 1-10:1

source ~/.bash_profile
# SGE_TASK_ID=1
rho=(0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5)
this_rho=${rho[$SGE_TASK_ID-1]}

simulation='Simulation_6_rhoCAV='$this_rho
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_rhoCAV_'$SGE_TASK_ID'_results'

for i in {1..1000};do
	python3 ../simulate.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt -n 1
	python3 ../caviar.py -l LD/ld_file_0.txt -z Z/z_file_0.txt -o $outfile'_'$i -r $this_rho
	rm -r LD
	rm -r Z

	speculate_file=$outfile'_'$i'_set.txt'
	python3 ../capture.py -s $speculate_file -t true_causal_set.txt -p '../'$simulation
done

mv $simulation'_config_size.txt' Config_Size
mv $simulation'_recall_rate.txt' Recall_Rate
