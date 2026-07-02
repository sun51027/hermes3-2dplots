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
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump 260421_reproduce_1e20_neumann_nowallpump -in 1e21 1e20 -o 260429_1e20_1e21_compare_benchmark -r omp -p outer_lower -m  benchmark
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump -in 1e21 -o 260429_1e21_benchmark -r omp -p outer_lower -m benchmark

############### UPDATED PLTOOL ################

# python3 make_plot.py -i 260217-cdn-46895-nowallpump_limiterOff_1e21 -in Feb_1e21 -o oldbranch_plot -r omp -p outer_lower -m special
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 -in Mar_1e21 -o oldbranch_plot -r omp -p outer_lower -m special
# python3 make_plot.py -i 260217-cdn-46895-nowallpump_limiterOff_1e21 -in Feb_1e21 -o 260506_balance -r omp -p outer_lower -m special
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 -in Mar_1e21 -o 260506_balance -r omp -p outer_lower -m special

# python3 make_plot.py -i 260114_cdn_46895_bndry_neumann  260217-cdn-46895-nowallpump_limiterOff_1e21 -in Jan_1.3e22 Feb_1e21 -o oldbranch_plot -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260421_reproduce_1e20_neumann_nowallpump 260417_reproduce_1e20_neumann_NV 260417_reproduce_1e20 -in neumann_nowallpump neumann decaylength -o 260430_bc_compare  -r omp -p outer_lower -m  multifiles
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump -o 260430_reproduce_1e20_neumann_nowallpump -in 1e21  -r omp -p outer_lower -m  benchmark  singlefile
# python3 make_plot.py -i 260430_reproduce_1e21_neumann_nowallpump_inoff  -o 260430_incoupling_OFF -in OFF    -r omp -p outer_lower -m  benchmark
# python3 make_plot.py -i 260430_reproduce_1e21_neumann_nowallpump_inoff 260427_reproduce_1e21_neumann_nowallpump -o 260430_incoupling -in OFF ON   -r omp -p outer_lower -m  benchmark
# python3 make_plot.py -i 260503_reproduce_5e20_neumann_nowallpump 260310_neutralimit_gcc12_9497476 -in reproduce_5e20 march_1e21 -o 260504_rep_mar_compare -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260421_reproduce_1e20_neumann_nowallpump 260503_reproduce_5e20_neumann_nowallpump 260427_reproduce_1e21_neumann_nowallpump -in 1e20 5e20 1e21 -o 260504_5e20_1e21_compare -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260503_reproduce_5e20_neumann_nowallpump 260504_reproduce_7e20_neumann_nowallpump -in 5e20 7e20 -o 260504_puff_compare -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260421_reproduce_1e20_neumann_nowallpump 260503_reproduce_5e20_neumann_nowallpump 260504_reproduce_7e20_neumann_nowallpump 260427_reproduce_1e21_neumann_nowallpump 260502_reproduce_5e21_neumann_nowallpump -in 1e20 5e20 7e20 1e21 5e21 -o 260504_puff_compare -r omp -p outer_lower -m special
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 -o 260506_balance -r omp -p outer_lower -m special
# python3 make_plot.py -i 260504_reproduce_7e20_neumann_nowallpump -o 260506_balance -r omp -p outer_lower -m special
################ Local run #####################

# python3 make_plot.py -i 260503_reproduce_5e20_neumann_nowallpump 260504_reproduce_7e20_neumann_nowallpump 260427_reproduce_1e21_neumann_nowallpump \
# -in 5e20 7e20 1e21 -o 260512_puff_compare -r omp -p outer_lower -m special
# python3 make_plot.py -i 260504_reproduce_7e20_neumann_nowallpump -o 260512_benchmark_7e20 -r omp -p outer_lower -m special
# python3 make_plot.py -i 260513_6e20_neumann_nowallpump_newgrid 260519_7e20_neumann_nowallpump 260519_1e21_neumann_nowallpump 260513_reproduce_6e20_neumann_nowallpump \
# -in 6e20 7e20 1e21 6e20_high -o 260520_lowresol_compr -r omp -p outer_lower -m benchmark special

# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260513_6e20_neumann_nowallpump_newgrid -in high low -o 260519_high_lowresol_6e20 -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260513_6e20_neumann_nowallpump_newgrid -o 260519_lowresol_6e20 -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260519_1e21_neumann_nowallpump 260513_reproduce_6e20_neumann_nowallpump -in Low_1e21 High_6e20 -o 260526_low_hig_comp -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260513_6e20_neumann_nowallpump_newgrid -o 260516_lowresol_6e20 -r omp -p outer_lower -m special
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump -o 260518_highresol_6e20 -r omp -p outer_lower -m special
# python3 make_plot.py -i 260503_reproduce_5e20_neumann_nowallpump -o 260518_highresol_6e20 -r omp -p outer_lower -m special
# python3 make_plot.py -i 260504_reproduce_7e20_neumann_nowallpump -o 260518_highresol_6e20 -r omp -p outer_lower -m special
# python3 make_plot.py -i 260518_re_6.5e20_from_6e20 260518_re_6.5e20_from_7e20 260513_reproduce_6e20_neumann_nowallpump -in 6.5e20_from_6e20 6.5e20_from_7e20 6e20 -o 260518_6.5e20_compare -r omp -p outer_lower -m special
# python3 make_plot.py -i 260503_reproduce_5e20_neumann_nowallpump 260513_reproduce_6e20_neumann_nowallpump 260504_reproduce_7e20_neumann_nowallpump -in 5e20 6e20 7e20 -o 260518_benchmark_all -r omp -p outer_lower -m special
# python3 make_plot.py -i 260503_reproduce_5e20_neumann_nowallpump 260513_reproduce_6e20_neumann_nowallpump 260504_reproduce_7e20_neumann_nowallpump -in 5e20 6e20 7e20 -o 260518_benchmark_all -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260526_5e20_neumann_nowallpump 260528_5e20_neumann_nowallpump_coreiz -in normal core_iz -o 260529_5e20_core_iz -r omp -p outer_lower -m special
# python3 make_plot.py -i 260526_5e20_neumann_nowallpump 260528_5e20_neumann_nowallpump_coreiz 260604_5e20_neumann_nowallpump_coreiz_scratch 260605_5e20_neumann_nowallpump_coreiz_a5 \
# -in normal core_iz from_scratch alpha_5 -o 260608_5e20_core_iz -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260526_5e20_neumann_nowallpump 260615_5e20_neumann_nowallpump_coreiz_izloss 260610_5e20_neumann_nowallpump_coreiz_onlyPF -in origin izloss_ON onlyPF -o 260615_coreiz -r omp -p outer_lower -m special benchmark
# python3 make_plot.py -i 260610_5e20_neumann_nowallpump_coreiz_NVd 260615_5e20_neumann_nowallpump_coreiz_izloss 260610_5e20_neumann_nowallpump_coreiz_onlyPF -in izloss_OFF izloss_ON onlyPF -o 260615_coreiz -r omp -p outer_lower -m special benchmark
# python3 make_plot.py -i 260601_reproduce_6e20_neumann_nowallpump_coreiz 260513_reproduce_6e20_neumann_nowallpump -in core_iz normal -o 260603_highresol_6e20_core_iz -r omp -p outer_lower -m special benchmark
# python3 make_plot.py -i 260617_5e20_corebc_return_True_izloss_True 260617_5e20_corebc_return_True_izloss_False 260617_5e20_corebc_izloss_True 260617_5e20_corebc_izloss_False 260526_5e20_neumann_nowallpump \
# -in return_O_izloss_O return_O_izloss_X return_X_izloss_O return_X_izloss_X normal -o 260617_coreiz_compr_all -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260519_1e21_neumann_nowallpump 260521_1e21_nowallpump_core_iz -in normal core_iz -o 270522_core_iz -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260619_lowrsl_5e20_scratch 260619_lowrsl_5e20_restart -in scratch restart -o 270621_5e20_scratch_restart -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260619_lowrsl_5e20_scratch 260619_lowrsl_5e20_corebc_OO_scratch 260619_lowrsl_5e20_corebc_XX_scratch \
# -in normal corebc_OO corebc_XX -o 270621_lowrsl_5e20_corebc_scratch -r omp -p outer_lower -m benchmark
# python3 make_plot.py -i 260619_lowrsl_5e20_restart 260619_lowrsl_5e20_corebc_OO_scratch 260619_lowrsl_5e20_corebc_XX_scratch \
#   -in normal corebc_OO corebc_XX -o 270621_lowrsl_5e20_corebc_restart -r omp -p outer_lower -m benchmark
python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260619_highrsl_6e20_corebc_return_X_izloss_X \
  -in normal corebc_XX -o 270701_highrsl_6e20_corebc_scratch -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260619_highrsl_6e20_corebc_return_O_izloss_O 260619_highrsl_6e20_corebc_return_X_izloss_X \
# -in normal corebc_OO corebc_XX -o 270621_highrsl_6e20_corebc_scratch -r omp -p outer_lower -m benchmark

######### Low and high resolution comparison
# python3 make_plot.py -i 260526_5e20_neumann_nowallpump 260513_reproduce_6e20_neumann_nowallpump -in Low_5e21 High_6e20 -o 260527_low_hig_comp -r omp -p outer_lower -m benchmark special
