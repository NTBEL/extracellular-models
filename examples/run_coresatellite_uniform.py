import numpy as np
import pandas as pd
from scipy.constants import N_A
from tifffile import imwrite as tiffwrite
from neuron.units import M
from extmodels.core_satellite import Model

model = Model(dx=10., xlo=-500., xhi=500., ylo=-500., yhi=500., zlo=-500., zhi=500.) 
model.radius_satellite = 20.
Qsst = 1.2e8
#r_inner = 0
r_core = model.radius_core
#inner_volume = (4/3)*np.pi*r_inner**3
core_volume = (4/3)*np.pi*r_core**3
initial_volume = core_volume 
volume_fraction = 0.2
sst_0 = ((Qsst / N_A) / initial_volume) / volume_fraction * 1e15 * M
model.sst.initial = lambda nd: sst_0 if (nd.x3d**2 + nd.y3d**2 + nd.z3d**2 < r_core**2) else 0
print('simulating...')
model.simulate(1, 50000, 200, nthread=16)
# Get and write out the z-projections as imagej tiff trajectories.
image2d_sensor = model.observables['zproject_mean_sensor_complex']
image2d_sst = model.observables['zproject_mean_sst']
dx = model.ecs._dx[0]
tiffwrite('trajectory_sensor.tif', np.array(image2d_sensor).astype(np.float32), imagej=True,
           metadata={'spacing' : dx, 'unit': 'micron', 'axes': 'TYX', 'fps':5})
tiffwrite('trajectory_sst.tif', np.array(image2d_sst).astype(np.float32), imagej=True,
           metadata={'spacing' : dx, 'unit': 'micron', 'axes': 'TYX', 'fps':5})
# Now get the core and satellite ROI mean values.
core_roi = model.observables['core_roi_mean']
satellite_roi = model.observables['satellite_roi_mean']
times = model.times
roi_data = pd.DataFrame({'time':times, 'core':core_roi, 'satellite':satellite_roi})
roi_data.to_csv('roi_data.csv')

print(len(model.times))

