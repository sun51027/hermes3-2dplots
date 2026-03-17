# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_pd0.nc --Nd 1.3e+22 --Pd 0.0 --puff omp --mode edit
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_v3.nc --Nd 1.3e+22 --puff omp --mode edit
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_corrected.nc --Nd 1.3e+22 --Pd 1.0e6 --puff omp --mode edit

#2026 new start 
# The original file Lloyd gave me is CDN_46895_v2.nc, but there was no Nd_src and Pd _src, so you need to assign value to it. You can't run it directly
# CDN_46895_corrected.nc was the modified one 
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_260108_1.3e21.nc --Nd 1.3e+21 --puff omp --mode edit
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_260105.nc --Nd 1.3e+22 --puff omp --mode edit
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_260120_nowallpump_2e21.nc --Nd 2.0e+21 --puff omp --mode edit
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_260129_nowallpump_1e21.nc --Nd 1.0e+21 --puff omp --mode edit
# python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/neutralrun/hermes-3/master/CDN_46895_260213_allwallpump_1e21.nc --Nd 1.0e+21 --puff omp --mode edit
python3 tune_puff.py --old-grid ~/neutralrun/hermes-3/master/CDN_46895_v2.nc --new-grid ~/CDN_46895_260314_nowallpump_3e20.nc --Nd 3.0e+20 --puff omp --mode edit
