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

simulation='Simulation_1_rho='$this_rho
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_rho_'$SGE_TASK_ID'_results'

for i in {1..1000};do
	python3 ../simulate_hardcode.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt
	python3 ../MCaviar.py -l LD -z Z -o $outfile'_multiple_'$i -r $this_rho
	python3 ../concatenate.py -z Z
	python3 ../caviar.py -l ../sample_data/50_LD.txt -z concatenated_z_score.txt -o $outfile'_original_'$i -r $this_rho
	rm -r LD
	rm -r Z

	speculate_file_m=$outfile'_multiple_'$i'_set.txt'
	speculate_file_o=$outfile'_original_'$i'_set.txt'
	python3 ../capture.py -s $speculate_file_m -t true_causal_set.txt -p '../'$simulation'_m'
	python3 ../capture.py -s $speculate_file_o -t true_causal_set.txt -p '../'$simulation'_o'
done

mv '../'$simulation'_m_config_size.txt' '../'Config_Size_m
mv '../'$simulation'_m_recall_rate.txt' '../'Recall_Rate_m
mv '../'$simulation'_o_config_size.txt' '../'Config_Size_o
mv '../'$simulation'_o_recall_rate.txt' '../'Recall_Rate_o