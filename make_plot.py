from pltool.plotting_single import run_single_plots
from pltool.plotting_multi import run_multi_plots 
import time, os


if __name__ == "__main__":
    
    start_time = time.time()
    run_single_plots()
    run_multi_plots()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    
