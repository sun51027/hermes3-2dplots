from pltool.plotting_single import * 
from pltool.plotting_multi_rings import * 
from pltool.plotting_multi_files import * 
from pltool.plotting_polygon import *
from pltool.benchmark import *
from pltool.common import *
import time, os


if __name__ == "__main__":
    
    start_time = time.time()

    setup_matplotlib()
    parser = build_base_parser()
    args = parser.parse_args()

    # Multiple files argument
    if len(args.input) > 1 and not args.input_name:
        parser.error("Please give short names for legend on multiple files comparison. \n")

    elif args.input_name and len(args.input) != len(args.input_name):
        parser.error(f"Number of file name is not correct, with {len(args.input)} files != {len(args.input_name)} name. \n")

    else:
        # If everything is correct, start to read Hermes-3
        max_len = max(len(name) for name in args.input)
        if max_len > 1 and args.input_name:
            print("Multiplfe files and its corresponding name:")
            for idx, f in enumerate(args.input):
                print(f"{idx}: {args.input[idx]:<{max_len}} -> {args.input_name[idx]}")
            print("\n")
            case = read_files(args.input, args.input_name)
        else:
            case = read_files(args.input)

    # # Read Hermes-3 case 
    # case = read_files(args.input, args.input_name)

    # Read SOLPS case 
    balance = read_solps_balance()

    ## Create output folder
    figures_png_path = args.output
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)


    ## Run program
    if "benchmark" in args.mode:
        print("========================================\n")
        print("Benchmark against single file mode")
        print("========================================")
        sepadd_array = [1]
        plot_bm_profiles_radial(case, balance, args.region_rad, figures_png_path)
        plot_bm_profiles_fieldline(case, balance, args.region_pol, sepadd_array, figures_png_path)

    if "multifiles" in args.mode:
        print("========================================\n")
        print("Compare multiple files mode")
        print("========================================")
        sepadd_array = 1
        plot_multifiles_profiles_fieldline(case, args.region_pol, sepadd_array, figures_png_path)
        plot_multifiles_profiles_radial(case, args.region_rad, figures_png_path)

    if "polygon" in args.mode:
        print("========================================\n")
        print("Plot two-dim polygon mode")
        print("========================================")

    if "singlefile" in args.mode:
        print("========================================\n")
        print("Plot profiles for single file (may contain several rings)")
        print("========================================")
        sepadd_array = [0, 1, 2, 3]
        plot_multi_profiles_fieldline(case, args.region_pol, sepadd_array, figures_png_path)
        plot_multi_profiles_radial(case, args.region_rad, figures_png_path)
        plot_plasma_overlap(case, args.region_pol, args.region_rad, figures_png_path)

    if "special" in args.mode:
        compare_last_timestep(case, balance,figures_png_path)
        plot_particle_balance_time(case, figures_png_path)

    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    
