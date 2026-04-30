#!/bin/bash
############################## OLD GRID - run in 10 cores

# python3 make_slurm_script.py -n 2 -c 5 -t 01:00:00 -d 26:03:20 -i /users/jpm590/scratch/260320_neutralimit_oldgrid_1e21 -j 260320_neutralimit_oldgrid_1e21 --restart
# python3 make_slurm_script.py -n 2 -c 5 -t 01:00:00 -d 26:03:20 -i /users/jpm590/scratch/260320_neutralimit_oldgrid_5e19 -j 260320_neutralimit_oldgrid_5e19  
# python3 make_slurm_script.py -n 2 -c 5 -t 01:00:00 -d 26:03:20 -i /users/jpm590/scratch/260322_neutralimit_oldgrid_1e21_petsc -j 260322_neutralimit_oldgrid_1e21_petsc --restart
# python3 make_slurm_script.py -n 2 -c 5 -t 01:00:00 -d 26:03:20 -i /users/jpm590/scratch/260322_neutralimit_oldgrid_1e21_petsc -j 260322_neutralimit_oldgrid_1e21_petsc --restart
############################## new GRID - run in 20 cores

#python3 make_slurm_script.py -n 1 -task 10 -c 2 -d 2025-09-30 -i /users/jpm590/2dspace/run/250929-shorter-test/ -j 2025-09-30-2D-MASU-restart --restart
#python3 make_slurm_script.py -n 2 -c 10 -d 2025-10-07 -i /users/jpm590/2dspace/run/251007-2D-MASTU/ -j 2025-10-07-2D-MASU-initial  -t 40:00:00
#python3 make_slurm_script.py -n 2 -c 10 -d 2025-10-07 -i /users/jpm590/2dspace/run/250929-shorter-test/ -j 2025-10-07-2D-MASU-restart --restart -t 40:00:00
#python3 make_slurm_script.py -n 2 -c 10 -d 2025-10-13 -i /users/jpm590/2dspace/run/251013-power-0.1MW/ -j 2025-10-13-power-0.1MW-initial  -t 40:00:00

#python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-1e21 -d 25:11:07 -j 2025-11-07-puff-1e21-initial -t 40:00:00 
#python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-1e22 -d 25:11:07 -j 2025-11-07-puff-1e22-initial -t 40:00:00 
#python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-5e20 -d 25:11:07 -j 2025-11-07-puff-5e20-initial -t 40:00:00 
#python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-5e22 -d 25:11:07 -j 2025-11-07-puff-5e22-initial -t 40:00:00 
# python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-1e21 -d 25:11:07 -j 2025-11-07-puff-1e21-restart -t 40:00:00 --restart 
# python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-1e22 -d 25:11:07 -j 2025-11-07-puff-1e22-restart -t 40:00:00 --restart 
# python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-5e20 -d 25:11:07 -j 2025-11-07-puff-5e20-restart -t 40:00:00 --restart 
# python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251107-tuned-puff-5e22 -d 25:11:07 -j 2025-11-07-puff-5e22-restart -t 40:00:00 --restart 

#python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/scratch/251112-tuned-puff-1e22 -d 25:11:12 -j 2025-11-12-puff-1e22-restart -t 40:00:00 --restart 
# python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/scratch/251117-pump-0.95-fix-Td -d 25:11:17 -j 2025-11-17-pump-fix-Td-restart -t 40:00:00 --restart 
#python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/2dspace/run/251112-tuned-puff-1e21 -d 25:11:12 -j 2025-11-12-puff-1e21-restart -t 40:00:00 --restart 
# python3 make_slurm_script.py -n 2 -c 10 -i /users/jpm590/scratch/251123-MASTU-cx-multiplier1000 -d 25:11:23 -j 2025-11-23-cx-multiplier1000 -t 40:00:00 --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 03:00:00 -d 26:03:19 -i /users/jpm590/scratch/260319_neutralimit_advec_compr_form -j 260319_pressure_form --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:03:19 -i /users/jpm590/scratch/260319_neutralimit_9497476_5e19 -j 260319_5e19 --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 00:10:00 -d 26:03:24 -i /users/jpm590/scratch/260324_neutralimit_flooring -j 260324_flooring --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:03:26 -i /users/jpm590/scratch/260325_neutralimit_cond_coeff -j 260326_cond_coeff --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:03:26 -i /users/jpm590/scratch/260326_neutral_test -j 260326_test --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:03:26 -i /users/jpm590/scratch/260327_neutral_bc -j 260327_bc --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:03:26 -i /users/jpm590/scratch/260327_neutral_flux0.1  -j 260327_flux0.1 --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 02:00:00 -d 26:03:30 -i /users/jpm590/scratch/260330_cx_floor -j 260330_cx_floor --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:04:01 -i /users/jpm590/scratch/260401_add_neutral_diagnose -j 260401_diagnose --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:04:01 -i /users/jpm590/scratch/260401_asym_puff -j 260401_asym_puff --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:04:02 -i /users/jpm590/scratch/260402_expand_puff -j 260402_expand_puff --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 02:00:00 -d 26:03:29 -i /users/jpm590/scratch/260327_neutral_lmax_0.01 -j 260329_lmax_0.01 --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 04:00:00 -d 26:03:27 -i /users/jpm590/scratch/260327_neutral_lmax  -j 260327_lmax --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:03:23 -i /users/jpm590/scratch/260323_neutralimit_decayBC_1e21 -j 260323_decayBC_1e21 --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:00:00 -d 26:04:07 -i /users/jpm590/scratch/260407_pump_95 -j 260407_pump_95 --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 01:30:00 -d 26:04:07 -i /users/jpm590/scratch/260407_fromscratch_1e21 -j 260407_fromscratch_1e21 
# python3 make_slurm_script.py -n 4 -c 5 -t 05:00:00 -d 26:04:10 -i /users/jpm590/scratch/260410_puff_1e19 -j 260410_puff_1e19 --restart  --mkdir
# python3 make_slurm_script.py -n 4 -c 5 -t 03:00:00 -d 26:04:08 -i /users/jpm590/scratch/260408_pump_95_originalBC -j  260408_pump_95_originalBC --restart 

###### updated make_slurm_script

# python3 make_slurm_script.py -n 4 -c 5 -t 03:00:00 -d 26:04:08 -i 260408_pump_95_originalBC -j  260408_pump_95_originalBC --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 03:00:00 -d 26:04:09 -i 260409_ion_coupling_off -j  260409_ion_coupling_off --restart 
# python3 make_slurm_script.py -n 4 -c 5 -t 05:00:00 -d 26:04:10 -i 260410_puff_1e19 -j 260410_puff_1e19 --restart  --mkdir
# python3 make_slurm_script.py -n 4 -c 5 -t 05:00:00 -d 26:04:10 -i 260410_puff_1e19 -j 260410_puff_1e19 --restart  --mkdir
# python3 make_slurm_script.py -n 4 -c 5 -t 02:00:00 -d 26:04:10 -i 260410_reproduce -j 260410_reproduce --restart  --mkdir
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:12 -i 260412_reproduce_1e19 -j 260412_reproduce_1e19 --restart  
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:17 -i 260417_reproduce_1e19 -j 260417_reproduce_1e19 --restart  
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:17 -i 260417_reproduce_1e20 -j 260417_reproduce_1e20 --restart  
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:17 -i 260417_reproduce_1e19_neumann -j 260417_reproduce_1e19_neumann --restart  
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:17 -i 260417_reproduce_1e19_neumann_NV -j 260417_reproduce_1e19_neumann_NV --restart  
# python3 make_slurm_script.py -n 2 -c 10 -t 10:00:00 -d 26:04:27 -i 260427_reproduce_1e21_neumann_nowallpump -j 260427_reproduce_1e21_neumann_nowallpump  --restart  
# python3 make_slurm_script.py -n 2 -c 10 -t 08:00:00 -d 26:04:21 -i 260421_reproduce_1e20_neumann_nowallpump -j 260421_reproduce_1e20_neumann_nowallpump  --restart  
# python3 make_slurm_script.py -n 2 -c 10 -t 11:00:00 -d 26:04:19 -i 260419_reproduce_1e20_neumann_NV_nowallpump -j 260419_reproduce_1e20_neumann_NV_nowallpump --restart  
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:17 -i 260417_reproduce_1e20_neumann_NV -j 260417_reproduce_1e20_neumann_NV --restart  
# python3 make_slurm_script.py -n 4 -c 5 -t 06:00:00 -d 26:04:17 -i 260417_reproduce_1e19_nowallpump -j 260417_reproduce_1e19_nowallpump --restart  
# python3 make_slurm_script.py -n 2 -c 10 -t 00:30:00 -d 26:04:29 -i 260429_reproduce_1e21_neumann_nowallpump_petsc -j 260429_reproduce_1e21_neumann_nowallpump_petsc --restart  
python3 make_slurm_script.py -n 2 -c 10 -t 02:00:00 -d 26:04:30 -i 260430_reproduce_1e21_neumann_nowallpump_inoff -j 260430_reproduce_1e21_neumann_nowallpump_inoff --restart  

#### TEST
# python3 make_slurm_script.py -n 1 -c 10 -t 00:10:00 -d 26:04:27 -i diffusion -j diffusion_test  
# python3 make_slurm_script.py -n 1 -c 20 -t 00:10:00 -d 26:04:27 -i 260427_slurm_test -j slurm_test 
