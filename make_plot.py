# from pltool.plotting_single import run_single_plots
from pltool.plotting_multi import run_multi_plots 
from pltool.plotting_multi_files import run_multifiles_plots 
from pltool.plotting_polygon import run_polygon 
from pltool.benchmark import *
from pltool.common import *
import time, os


if __name__ == "__main__":
    
    start_time = time.time()

    setup_matplotlib()
    parser = build_base_parser()
    args = parser.parse_args()
    ## Read Hermes-3 case 
    case = read_files(args.input)
    ## Read SOLPS case 
    balance = read_solps_balance()

    ## output folder
    figures_png_path = args.output
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)

    if "benchmark" in args.mode:
        print("Benchmark against single file mode")
        print("========================================")
        sepadd_array = [1]
        plot_bm_profiles_radial(case, balance, args.region_rad, figures_png_path)
        plot_bm_profiles_fieldline(case, balance, args.region_pol, sepadd_array, figures_png_path)

    if "multifiles" in args.mode:
        print("Compare multiple files mode")
        print("========================================")
    if "polygon" in args.mode:
        print("Plot two-dim polygon mode")
        print("========================================")
    if "multi_benchmark" in args.mode:
        print("Benchmark against multiple files mode")
        print("========================================")
    # run_multi_plots()
    # run_multifiles_plots()
    # run_polygon()
    # run_benchmark()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    
