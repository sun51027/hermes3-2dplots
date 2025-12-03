from pltool.plotting_single import main #plot_single_profiles
import time, os


if __name__ == "__main__":
    
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    
