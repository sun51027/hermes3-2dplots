import boutdata
from boutdata import restart
from boutdata import collect

original_path = "/users/jpm590/scratch/260213-cdn-46895-wallRecyON_1e21/"
new_path = "/users/jpm590/scratch/260217-cdn-46895-wallRecyON_1e21/"
restart.redistribute(20,  nxpe=1, path=original_path, output=new_path)

# Old restart file
# print("Old restart file:")
# nd = collect("Nd+", path="/users/jpm590/scratch/260207-cdn-46895-nowallpump_1e21/", info=True)
# print(nd.shape)

# print("new restart file:")
# nd = collect("Nd+", path="/users/jpm590/scratch/260217-cdn-46895-nowallpump_1e21/",prefix="BOUT.restart")
# print(nd.shape)
