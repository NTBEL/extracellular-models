import numpy as np
import pandas as pd
from scipy.constants import N_A
from tifffile import imwrite as tiffwrite
from neuron.units import M, nM
from extmodels.dynorphin_pointsource import model

binding_site_conc = np.logspace(0, 4, 5, endpoint=True)
binding_kd = binding_site_conc.copy()
print(list(binding_site_conc))
print(list(binding_kd))
dx = model.ecs._dx[0]
for bs_conc in binding_site_conc:
    model.sensor.initial = bs_conc * nM
    for b_kd in binding_kd:
        model.sensor_kd = b_kd * nM
        print("bs_conc: ", bs_conc, " b_kd: ", b_kd, " simulating...")
        model.simulate(5, 4000, 20, nthread=4)
        # Get and write out the z-projections as imagej tiff trajectories.
        image2d_dyn = model.observables["zproject_mean_dyn"]
        tiffwrite(
            "./image_trajs/trajectory_zproject_dynorphin_bsc{}_bkd{}.tif".format(
                bs_conc, b_kd
            ),
            np.array(image2d_dyn).astype(np.float32),
            imagej=True,
            metadata={"spacing": dx, "unit": "micron", "axes": "TYX", "fps": 10},
        )
        image2d_complex = model.observables["zproject_mean_sensor_complex"]
        tiffwrite(
            "./image_trajs/trajectory_zproject_sensor-complex_bsc{}_bkd{}.tif".format(
                bs_conc, b_kd
            ),
            np.array(image2d_complex).astype(np.float32),
            imagej=True,
            metadata={"spacing": dx, "unit": "micron", "axes": "TYX", "fps": 10},
        )
        quit()
