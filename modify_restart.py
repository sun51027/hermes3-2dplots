import boutdata
from boutdata import restart
from boutdata import collect


def change_n_core():
    original_path = "/users/jpm590/scratch/260213-cdn-46895-wallRecyON_1e21/"
    new_path = "/users/jpm590/scratch/260217-cdn-46895-wallRecyON_1e21/"
    restart.redistribute(20,  nxpe=1, path=original_path, output=new_path)

def change_gridfile():
    old_grid = "/users/jpm590/gridfile/CDN_46895_260129_nowallpump_1e21.nc"
    new_grid = "/users/jpm590/gridfile/mu1af6-tunepuff.nc"
    restartfile = "/users/jpm590/260320_neutralimit_oldgrid_1e21/"
    output_path = "/users/jpm590/260320_neutralimit_oldgrid_1e21/"
    restart.change_grid( old_grid, new_grid, restartfile, output_path)


change_gridfile()
# Old restart file
# print("Old restart file:")
# nd = collect("Nd+", path="/users/jpm590/scratch/260207-cdn-46895-nowallpump_1e21/", info=True)
# print(nd.shape)

# print("new restart file:")
# nd = collect("Nd+", path="/users/jpm590/scratch/260217-cdn-46895-nowallpump_1e21/",prefix="BOUT.restart")
# print(nd.shape)
