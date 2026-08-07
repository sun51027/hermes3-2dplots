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
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260619_highrsl_6e20_corebc_return_X_izloss_X \
#   -in normal corebc_XX -o 270701_highrsl_6e20_corebc_scratch -r omp -p outer_lower_sol -m benchmark

################ Updated sdtool #####################
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260619_highrsl_6e20_corebc_return_O_izloss_O 260619_highrsl_6e20_corebc_return_X_izloss_X \
# -in normal corebc_OO corebc_XX -o 260703_corebc_OO_XX -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260504_reproduce_7e20_neumann_nowallpump 260705_highrsl_7e20_corebc_XX 260705_highrsl_7e20_corebc_OO \
#   -in normal corebc_XX corebc_OO -o 260717_7e20_corebc_OO_XX -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260504_reproduce_7e20_neumann_nowallpum 260705_highrsl_7e20_corebc_XX 260705_highrsl_7e20_corebc_OO \
#   -in normal_6e20 normal_7e20 core_XX_7e20 corebc_OO_7e20 -o 260706_7e20_corebc_OO_XX -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260703_highrsl_6e20_corebc_XX_alpha_10 260619_highrsl_6e20_corebc_return_X_izloss_X \
# -in normal XX_alpha_10 XX_1 -o 260703_nalpha -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260427_reproduce_1e21_neumann_nowallpump 260619_highrsl_6e20_corebc_return_X_izloss_X \
#   -in 6e20 1e21 6e20_coreXX -o 260706_6e20_1e21 -r omp -p outer_lower_sol -m benchmark
#
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260706_highrsl_1e21 260706_highrsl_1e21_corebc_OO 260706_highrsl_1e21_corebc_XX \
# -in 6e20 1e21 1e20_coreOO 1e21_coreXX -o 260707_6e20_1e21_corebc -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump 260706_highrsl_1e21 260706_highrsl_1e21_corebc_OO \
# -in 6e20 1e21 1e21_coreOO -o 260707_6e20_1e21_corebc -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260619_highrsl_6e20_corebc_return_X_izloss_X 260706_highrsl_1e21 260706_highrsl_1e21_corebc_OO \
#   -in 6e20_XX 1e21 1e21_coreOO -o 260707_6e20_1e21_corebc_v2 -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260706_highrsl_1e21 260707_lowrsl_1e21 \
#   -in high_1e21 Low_1e21 -o 260716_rsl_compr_1e21 -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260707_lowrsl_1e21_core_XX 260707_lowrsl_1e21_core_OO 260707_lowrsl_1e21 \
#   -in 1e21_XX 1e21_OO 1e21_normal -o 260715_lowrsl_1e21_corebc -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260706_highrsl_1e21 260706_highrsl_1e21_corebc_OO 260706_highrsl_1e21_corebc_XX \
#   -in 1e21_normal 1e21_OO 1e21_XX -o 260715_highrsl_1e21_corebc_lucinai -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260619_lowrsl_5e20_restart 260513_reproduce_6e20_neumann_nowallpump -in 5e20_low 6e20_high -o 260717_low5e20_high6e20 -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_cx10 -in core_iz core_iz_cx10 -o 260721_lowrsl_6e20_cx10 -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_nlmx100 \
#   -in core_iz core_iz_nlmx100 \
#   -o 260805_lowrsl_6e20_nlmx100 \
#   -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20 -o 260723_lowrsl_6e20_forTAP -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20 260717_lowrsl_6e20_coreOO -in normal core_iz -o 260723_lowrsl_6e20_forTAP -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_highrsl_6e20 260706_highrsl_1e21 -in 6e20 1e21 -o 260806_highrsl_6e20_1e21_origin -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_cx100 -in 6e20_coreiz 6e20_coreiz_cx100 -o 260723_lowrsl_6e20_solpscx100 -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_cx10 260717_lowrsl_6e20_coreOO_cx0.1 -in normal cx10 cx0.1 -o 260724_lowrsl_6e20cx -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20 260707_lowrsl_1e21 -in 6e20 1e21 -o 260723_lowrsl_6e20_1e21_forTAP  -m benchmark
# python3 make_plot.py -i 260717_highrsl_6e20 260706_highrsl_1e21 -in 6e20 1e21 \
#  -si SOLPS_Luciani_off -sin 1e21 \
#  -o 260729_highrsl_6e20_1e21_forPR  -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_cx0.1 260717_lowrsl_6e20_coreOO_cx100 -in normal coreiz coreiz_cx0.1 coreiz_cx100 \
# -o 260728_lowrsl_6e20_cx -r omp -p outer_lower_sol -m benchmark
# python3 make_plot.py -i 260717_lowrsl_6e20 260717_lowrsl_6e20_coreOO \
#  -in normal coreiz  \
#   -o 260803_lowrsl_6e20_coreiz_polygon \
#   -m benchmark \
#   --solps
# python3 make_plot.py -i 260717_highrsl_6e20 260717_highrsl_6e20_coreOO \
#  -in normal coreiz  \
#   -o 260803_highrsl_6e20_coreiz_rings \
#   -m benchmark \
#   --solps
# python3 make_plot.py -i 260717_lowrsl_6e20 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_cx0.1 260717_lowrsl_6e20_coreOO_cx100 -in normal coreiz coreiz_cx0.1 coreiz_cx100 \
# -o 260728_lowrsl_6e20_solps_cx100_1e22 -r omp -p outer_lower_sol -m benchmark
python3 make_plot.py -i 260806_highrsl_6e20_solps_1e21_plasmainit -in fix_plasma_bg -o 260807_fix_pls_bg -m benchmark

######### Multiple SOLPS files ########
# python3 make_plot.py -m benchmark \
#   -i 260717_lowrsl_6e20_coreOO 260717_lowrsl_6e20_coreOO_cx100 \
#   -in 6e20 6e20_cx100\
#   -si SOLPS_Luciani_off_cx100 SOLPS_Luciani_off_cx100_1e22 \
#   -sin 1e21_cx100 1e22_cx100 \
#   -o 260729_lowrsl_6e20_solps_cx100
#260717_lowrsl_6e20_coreOO_cx0.1 260717_lowrsl_6e20_coreOO_cx100 \
#coreiz_cx0.1 coreiz_cx100 \

############## Polygon ###############
# python3 make_plot.py -i 260513_reproduce_6e20_neumann_nowallpump -o 260703_solps_polygon -m polygon --solps
# python3 make_plot.py -i 260707_lowrsl_1e21 -o 260706_1e21_polygon -m polygon --solps
# python3 make_plot.py -i 260427_reproduce_1e21_neumann_nowallpump -o 260706_1e21_polygon -m polygon --solps
