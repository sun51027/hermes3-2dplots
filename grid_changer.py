import xhermes
from xhermes import *

hp = HypnotoadGrid(gridfilepath="../hypnotoad/CDN_46895_0511.nc")
# hp.check_decomposition(nxpe=10, nprocs=3)
# Call the function and save the returned results
is_valid, message = hp.check_decomposition(nxpe=1, nprocs=10)
# Print the results so you can see them
print(f"Is valid: {is_valid}")
print(f"Message: {message}")
