import numpy as np
import pandas as pd
from scipy.constants import N_A
from tifffile import imwrite as tiffwrite
from neuron.units import M
from extmodels.calcein_pointsource_asymm import Model

tort_vals_neg = np.linspace(1.,2.5, 16, endpoint=True)
tort_vals_pos = tort_vals_neg.copy()
print(list(tort_vals_neg))
print(list(tort_vals_pos))
steady_vals = list()
for tneg in tort_vals_neg:
    for tpos in tort_vals_pos:
        model = Model(tortuosity_neg=tneg, tortuosity_pos=tpos) 
        print('neg: ',tneg, 'pos: ', tpos, ' simulating...')
        model.simulate(5, 4000, 20, nthread=8)
        # Get and write out the z-projections as imagej tiff trajectories.
        image2d = model.observables['zproject_mean_calcein']
        dx = model.ecs._dx[0]
        tiffwrite('./image_trajs/trajectory_zproject_calcein_n{}_p{}.tif'.format(tneg, tpos), np.array(image2d).astype(np.float32), imagej=True,
                   metadata={'spacing' : dx, 'unit': 'micron', 'axes': 'TYX', 'fps':10})
        # Now get the core and satellite ROI mean values.
        num_neg = model.observables['num_neg'][1:]
        num_pos = model.observables['num_pos'][1:]
        ratio = num_pos / num_neg
        times = model.times[1:]
        out_data = pd.DataFrame({'time':times, 'Nneg':num_neg, 'Npos':num_pos, 'ratio':ratio})
        out_data.to_csv('./out_csv/out_data_n{}_p{}.csv'.format(tneg, tpos))
        print("    steady-state ratio (pos/neg): {}".format(np.round(ratio[-1], 1)))
        steady_vals.append([tneg, tpos, ratio[-1]])
steady_vals = np.array(steady_vals)
steady_data = pd.DataFrame({'tort_neg':steady_vals[::,0], 'tort_pos':steady_vals[::,1], 'ratio':steady_vals[::,2]})
steady_data.to_csv('steady_state_ratios.csv') 
