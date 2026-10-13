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
# python3 make_plot.py -i 260717_lowrsl_6e20 260707_lowrsl_1e21 -in 6e20 1e21 -o 260723_lowrsl_6e20_1e21_forTAP -m benchmark
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
# python3 make_plot.py \
#   -i 260717_highrsl_6e20 260717_highrsl_6e20_coreOO \
#   -in original ionising_core \
#   -si SOLPS_origin_1e21 \
#   -sin 1e21 \
#   -o 261005_benchmark_poster -m benchmark
# python3 make_plot.py -i 260310_neutralimit_gcc12_9497476 -in restart_dep -sin . -o 260824_restart_dep_forPR -m polygon --solps
# python3 make_plot.py -i 260806_highrsl_6e20_solps_1e21_plasmainit -in fix_plasma_bg -o 260907_fix_pls_bg_Spol -m benchmark
# python3 make_plot.py -i 260916_highrsl_fixplsbg_fromsolps -in fix_plasma_bg -o 260918_fix_pls_bg_Spar -m benchmark
# python3 make_plot.py -i 260916_highrsl_fixplsbg_fromsolps_cx100_1e21 -in fix_plasma_bg -si SOLPS_Luciani_off_cx100 -sin SOLPS_cx100_1e21 -o 260917_fix_pls_bg_Spar_cx_1e21 -m benchmark
# python3 make_plot.py -i 260916_highrsl_fixplsbg_fromsolps_cx100_1e22 -in fix_plasma_bg -si SOLPS_Luciani_off_cx100_1e22 -sin SOLPS_cx100_1e22 -o 260917_fix_pls_bg_Spar_cx_1e22 -m benchmark
# python3 make_plot.py -i 260916_highrsl_fixplsbg_fromsolps 260916_highrsl_fixplsbg_fromsolps_cx100_1e21 260916_highrsl_fixplsbg_fromsolps_cx100_1e22 \
#     -in fixbg fixbg_cx100_1e21 fixbg_cx100_1e22  \
#     -sin origin cx100_1e21 cx100_1e22 \
#     -si SOLPS_Luciani_off SOLPS_Luciani_off_cx100 SOLPS_Luciani_off_cx100_1e22 \
#     -o 260918_fix_pls_bg_all -m benchmark
# python3 make_plot.py -i 260916_highrsl_fixplsbg_fromsolps 260717_highrsl_6e20_coreOO \
#   -in fixbg_1e21 6e20_coreiz \
#   -sin origin_1e21 \
#   -si SOLPS_Luciani_off \
#   -o 260918_fix_pls_bg_vs_6e20_coreiz -m benchmark
# python3 make_plot.py -i 260916_highrsl_fixplsbg_fromsolps 260916_highrsl_fixplsbg_fromsolps_cx100_1e21 \
#   -in fixbg_1e21 fixbg_cx100_1e21 \
#   -si SOLPS_Luciani_off SOLPS_Luciani_off_cx100 \
#   -sin origin_1e21 cx100_1e21 \
#   -o 260918_fix_pls_bg_vs_origin_cx100_1e21 -m benchmark
# python3 make_plot.py -i 260921_highrsl_fixplsbg_fromsolps_1e21 260707_lowrsl_1e21_core_OO \
#   -in fixbg_1e21 lowrsl_origin_1e21\
#   -si SOLPS_Luciani_off \
#   -sin origin_1e21 \
#   -o 260921_fixplsbg_1e21_vs_lowrsl_1e21 -m benchmark
# python3 make_plot.py -i 260921_highrsl_fixplsbg_fromsolps_1e21 260717_highrsl_6e20_coreOO \
#   -in fixbg_1e21 origin_6e20\
#   -si SOLPS_Luciani_off \
#   -sin origin_1e21 \
#   -o 260921_fixplsbg_1e21_vs_normal_6e20 -m benchmark
# python3 make_plot.py -i 260921_highrsl_fixplsbg_fromsolps_1e21 \
#   -in fixbg_1e21 \
#   -si SOLPS_Luciani_off \
#   -sin origin_1e21 \
#   -o 260921_fixplsbg_vs_1e21_coreiz -m benchmark


# python3 make_plot.py -i 260921_highrsl_fixplsbg_fromsolps_1e21 260921_highrsl_fixplsbg_fromsolps_cx100_1e22_hermes_cx100 \
#   -in fixbg_1e21 fixbg_cx100_1e22\
#   -si SOLPS_Luciani_off SOLPS_Luciani_off_cx100_1e22 \
#   -sin origin_1e21 cx100_1e22 \
#   -o 260921_fixplsbg_1e21_vs_cx100_1e22 -m benchmark
# python3 make_plot.py -i 260921_highrsl_fixplsbg_fromsolps_1e21 260923_highrsl_fixplsbg_fromsolps_1e21_roomT \
#   -in fixbg_1e21 fixbg_1e21_roomT \
#   -si SOLPS_Luciani_off \
#   -sin origin_1e21 \
#   -o 260923_fixplsbg_1e21_roomT -m benchmark

# python3 make_plot.py -i 260923_highrsl_fixplsbg_fromsolps_1e21_roomT 260923_highrsl_fixplsbg_fromsolps_1e21_roomT_nlmx0.01 260923_highrsl_fixplsbg_fromsolps_1e21_roomT_nlmx100 \
#   -in nlmx0.1 nlmx0.01 nlmx100\
#   -si SOLPS_origin_1e21 \
#   -sin origin_1e21 \
#   -o 260924_fixplsbg_1e21_nlmx -m benchmark

# python3 make_plot.py -i 260923_highrsl_fixplsbg_fromsolps_1e21_roomT 261001_highrsl_fixplsbg_fromsolps_1e21_roomT_gd \
#   -in original fix_guardcell\
#   -si SOLPS_origin_1e21 \
#   -sin origin_1e21 \
#   -o 261001_fixplsbg_guardcell -m benchmark

# python3 make_plot.py -i 261001_highrsl_fixplsbg_fromsolps_1e21_roomT_gd \
#   -in fixed_plasma_1e21\
#   -si SOLPS_origin_1e21 \
#   -sin 1e21\
#   -o 261005_fixplsbg_poster\
#   -m benchmark
# python3 make_plot.py -i  261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd 260923_highrsl_fixplsbg_fromsolps_cx100_1e21_hermes_cx100_roomT \
#   -in 1e21_cx100_gd 1e21_cx100 \
#   -si SOLPS_1e21_cx100 \
#   -sin 1e21_cx100 \
#   -o 261013_1e21_cx100_gdfix\
#   -m benchmark
# python3 make_plot.py -i  261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd \
#   -in 1e21_cx100 \
#   -si SOLPS_1e21_cx100 \
#   -sin 1e21_cx100 \
#   -o 261013_1e21_cx100\
#   -m benchmark polygon \
#   --solps \
#   --same_limit 
# python3 make_plot.py -i 261001_highrsl_fixplsbg_fromsolps_cx100_1e22_roomT_gd \
#   -in 1e22_cx100 \
#   -si SOLPS_1e22_cx100 \
#   -sin 1e22_cx100 \
#   -o 261013_1e22cx100\
#   -m benchmark polygon \
#   --solps \
#   --same_limit
python3 make_plot.py -i 261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd 261001_highrsl_fixplsbg_fromsolps_cx100_1e22_roomT_gd \
  -in 1e21_cx100 1e22_cx100 \
  -si SOLPS_1e21_cx100 SOLPS_1e22_cx100 \
  -sin 1e21_cx100 1e22_cx100 \
  -o 261013_all_cx100\
  -m benchmark

# python3 make_plot.py -i 260706_highrsl_1e21_corebc_OO 260717_highrsl_6e20_coreOO \
#   -in 1e21 6e20 \
#   -si SOLPS_origin_1e21 \
#   -o 261013_highrsl_1e21_6e20_preview \
#   -m benchmark

# python3 make_plot.py -i 260717_highrsl_6e20_coreOO 260717_highrsl_6e20\
#   -in 6e20_coreOO 6e20 \
#   -si SOLPS_origin_1e21 \
#   -o 261013_highrsl_6e20_polygon \
#   -m polygon \
#   --solps \
#   --same_limit
# python analyse_penetration.py --case 261001_highrsl_fixplsbg_fromsolps_1e21_roomT_gd --solps ../SOLPS_origin_1e21 -o 261005_fixplsbg_poster
# python analyse_penetration.py --case 261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd --solps ../SOLPS_1e21_cx100 -o 261013_1e21_cx100
# python analyse_penetration.py --case 261001_highrsl_fixplsbg_fromsolps_cx100_1e22_roomT_gd --solps ../SOLPS_1e22_cx100 -o 261013_1e22_cx100

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

# python3 make_plot.py -i 260923_highrsl_fixplsbg_fromsolps_1e21_roomT \
#   -si SOLPS_Luciani_off \
#   -o 260923_fixplsbg_1e21_roomT \
#   -m polygon \
#   --solps \
#   --same_limit
