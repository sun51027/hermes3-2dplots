from pltool.plotting_single import run_single_plots
from pltool.plotting_multi import run_multi_plots 
from pltool.plotting_multi_files import run_multifiles_plots 
from pltool.plotting_polygon import run_polygon 
from pltool.benchmark import run_benchmark
import time, os


if __name__ == "__main__":
    
    start_time = time.time()
    # run_single_plots()
    # run_multi_plots()
    # run_multifiles_plots()
    # run_polygon()
    run_benchmark()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    
