# python3 single_plot.py -i 251007-2D-MASTU -o 251007_origin_omp_outerlower-ring15 -r omp -p outer_lower --sepadd 15
# python3 single_plot.py -i 251007-2D-MASTU -o 251007_origin_xtarget_outerlower-ring5 -r outer_lower_target -p outer_lower --sepadd 5
# python3 single_plot.py -i 251105-pump-0.95 -o 251105-pump-0.95_omp_outerlower-ring15 -r omp -p outer_lower --sepadd 15
# python3 single_plot.py -i 251105-pump-0.95 -o 251105-pump-0.95_xtarget_outerlower-ring5 -r outer_lower_target -p outer_lower --sepadd 5


# python3 single_plot.py -i 251105-pump-0.95 -o 251105-pump-0.95_omp_outerlower-log -r omp -p outer_lower --sepadd 1 -s log
# python3 single_plot.py -i 251105-pump-0.95 -o 251105-pump-0.95_omp_outerlower -r omp -p outer_lower --sepadd 1

############ Files moved to scratch #############

# python3 single_plot.py  -r omp -p outer_lower --sepadd 1 --scale log    -i 251007-2D-MASTU -o MASTU-origin-log
# python3 single_plot.py  -r omp -p outer_lower --sepadd 1 --scale linear -i 251007-2D-MASTU -o MASTU-origin-linear
# python3 single_plot.py -i 251112-tuned-puff-1e22 -o 251112-1e22 -r omp -p outer_lower --sepadd 1 
# python3 single_plot.py  -r omp -p outer_lower --sepadd 1 --scale log    -i 251119-MASTU-newbranch-rerun -o MASTU-neutralrun-log
# python3 single_plot.py -i 251119-MASTU-newbranch-rerun -o 251119-mastu-newbranch-log -r omp -p outer_lower --sepadd 1 --scale log
# python3 single_plot.py -i 251119-MASTU-newbranch-rerun -o 251119-mastu-newbranch-log -r omp -p outer_lower --sepadd 1 --scale log
# python3 multi_plot.py -i 251119-MASTU-newbranch-rerun -o 251119-mastu-newbranch-log -r omp -p outer_lower --sepadd 1 --scale log


# python3 make_plot.py -i 251212-cdn-46895-corrected -o 251212-newgrid-corrected -r omp -p outer_lower
# python3 make_plot.py -i 260105-cdn-46895-david-param -o 260105-cdn-46895-david-param -r omp -p outer_lower
# python3 make_plot.py -i 251205-cdn-46895-old-param -o 251212-newgrid-old-param -r omp -p outer_lower
# python3 multi_plot.py -i 251007-2D-MASTU -o 251201-mastu-original-multi -r omp -p outer_lower
# python3 multi_plot.py -i 251123-MASTU-cx-multiplier1000 -o 251201-mastu-cx-multiplier1000-multi -r omp -p outer_lower

# python3 make_plot.py -i 260112-cdn-46895-david-param -o 260123-cdn-46895-david-param -r omp -p outer_lower
# python3 make_plot.py -i 260114-cdn-46895-bndry-neumann -o 260123-cdn-46895-bndry-neumann -r omp -p outer_lower
# python3 make_plot.py -i 260207-cdn-46895-nowallpump_1e21 -o 260207-cdn-46895-nowallpump_1e21 -r omp -p outer_lower

# python3 make_plot.py -i 260217-cdn-46895-nowallpump_limiterOff_1e21 -o 260217-nowallpump_limiterOff_1e21 -r omp -p outer_lower
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 -o 260331_neutralimit_gcc12_9497476 -r omp -p outer_lower
# python3 make_plot.py -i 260314_neutralimit_9497476_3e20 -o 260314_9497476_cvode_3e20 -r omp -p outer_lower
# python3 make_plot.py -i 260316_neutralimit_9497476_3e20_from0217 -o 260316_9497476_cvode_3e20_from0217 -r omp -p outer_lower
# python3 make_plot.py -i 260217-cdn-46895-nowallpump_limiterOff_1e21_extend -o 260316-cdn-46895-nowallpump_limiterOff_1e21_extend -r omp -p outer_lower
# python3 make_plot.py -i 260319_neutralimit_advec_compr_form -o  260319_neutralimit_advec_compr_form  -r omp -p outer_lower
# python3 make_plot.py -i 260319_neutralimit_9497476_5e19  -o 260319_neutralimit_9497476_5e19 -r omp -p outer_lower
# python3 make_plot.py -i 260327_neutral_flux0.1 -o 260327_neutral_flux0.1 -r omp -p outer_lower
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476  -o 260310_neutralimit_gcc12_9497476 -r omp -p outer_lower
# python3 make_plot.py -i 260407_pump_95  -o 260410_pump_95 -r omp -p outer_lower 
# python3 make_plot.py -i 260409_ion_coupling_off  -o 260410_ion_coupling_off -r omp -p outer_lower 
# python3 make_plot.py -i 260410_puff_1e19  -o 260410_puff_1e19 -r omp -p outer_lower 


################# Compare different simulation ######################
# python3 make_plot.py -i 260207* -o 260207_multiple_comparison -r omp -p outer_lower 
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 260407_pump_95  -o 260413_pump95_compare -r omp -p outer_lower 
# python3 make_plot.py -i 260410_puff_1e19 260310_neutralimit_gcc12_9497476 -o 260410_puff_compare -r omp -p outer_lower 
# python3 make_plot.py -i 260217-cdn-46895-nowallpump_limiterOff_1e21 260310_neutralimit_gcc12_9497476 -o 260410_oldnew -r omp -p outer_lower 


################# Benchamrk against solps ######################
# python3 make_plot.py -i 260305-nowallpump-limiterOff-1e21-cvode  -o 260310-cvode  -r omp -p outer_lower
# python3 make_plot.py -i 260217-cdn-46895-nowallpump_limiterOff_1e21   -o 260407_new_bal  -r omp -p outer_lower
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 -o 260407_new_bal_neutralimit -r omp -p outer_lower
# python3 make_plot.py -i 260327_neutral_lmax -o 260328_neutral_lmax_benchmark -r omp -p outer_lower
# python3 make_plot.py -i 260417_reproduce_1e20 260417_reproduce_1e20_neumann_NV 260419_reproduce_1e20_neumann_NV_nowallpump -o 260421_reproduce0417_compare -r omp -p outer_lower

######## New SOLPS fil
# python3 make_plot.py -i 260407_pump_95 -o 260407_pump_95_new -r omp -p outer_lower
# python3 make_plot.py -i 260410_reproduce -o 260414_reproduce_new -r omp -p outer_lower
# python3 make_plot.py -i 260412_reproduce_1e19 -o 260417_reproduce_1e19_new -r omp -p outer_lower
# python3 make_plot.py -i 260417_reproduce_1e19 -o 260417_reproduce_1e19 -r omp -p outer_lower
# python3 make_plot.py -i 260417_reproduce_1e20 -o 260418_reproduce_1e20 -r omp -p outer_lower
# python3 make_plot.py -i 260417_reproduce_1e19 -o 260418_reproduce_1e19 -r omp -p outer_lower
# python3 make_plot.py -i 260419_reproduce_1e20_neumann_NV_nowallpump -o 260419_reproduce_1e20_neumann_NV_nowallpump -r omp -p outer_lower
# python3 make_plot.py -i 260417_reproduce_1e20_neumann_NV -o 260421_reproduce_1e20_neumann_NV -r omp -p outer_lower
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump -o 260427_reproduce_1e21_neumann_nowallpump -r omp -p outer_lower
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump 260421_reproduce_1e20_neumann_nowallpump -in 1e21 1e20 -o 260429_1e20_1e21_compare -r omp -p outer_lower -m  multifiles
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump 260421_reproduce_1e20_neumann_nowallpump -in 1e20 1e21 -o 260429_1e20_1e21_compare_benchmark -r omp -p outer_lower -m  benchmark
python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump -in 1e21 -o 260429_1e21_benchmark -r omp -p outer_lower -m benchmark 


# python3 make_plot.py -i 260407_pump_95 -o 260418_pump_95 -r omp -p outer_lower