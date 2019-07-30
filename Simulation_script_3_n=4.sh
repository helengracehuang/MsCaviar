#!/bin/bash
# $ -S /bin/bash
# $ -N job-simulation_mCAVIAR_rosemary_1
# $ -cwd = run from this current working directory (relative paths)
# -o stdout-Simulation_3_numOfSnps.out
# $ -l h_data=12G,h_rt=24:00:00
# $ -t 1-1000:1

source ~/.bash_profile
# SGE_TASK_ID=3
this_num_of_snps=4

sim='Simulation_3_numOfSnps='$this_num_of_snps
simulation='Simulation_3_numOfSnps='$this_num_of_snps'_'$SGE_TASK_ID
echo $simulation
mkdir $simulation
cd $simulation

outfile='sim_numOfSnps_4_'$SGE_TASK_ID'_results'

python3 ../simulate.py -l ../sample_data/50_LD.txt -z ../sample_data/50_Z.txt -c $this_num_of_snps
python3 ../MCaviar.py -l LD -z Z -o $outfile'_'$i -c $this_num_of_snps
# rm -r LD
# rm -r Z

speculate_file=$outfile'_'$i'_set.txt'
python3 ../capture.py -s $speculate_file -t true_causal_set.txt -p '../'$sim

# mv '../'$sim'_config_size.txt' '../'Config_Size
# mv '../'$sim'_recall_rate.txt' '../'Recall_Rate

mv '../'$simulation ../Simulation_3_numOfSnps=4
