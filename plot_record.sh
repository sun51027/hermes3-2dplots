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


python3 make_plot.py -i 251205-cdn-46895-old-param -o 251212-newgrid-old-param -r omp -p outer_lower
# python3 multi_plot.py -i 251007-2D-MASTU -o 251201-mastu-original-multi -r omp -p outer_lower
# python3 multi_plot.py -i 251123-MASTU-cx-multiplier1000 -o 251201-mastu-cx-multiplier1000-multi -r omp -p outer_lower
