# extracellular-models

This package contains a set of reaction-diffusion models to simulate the dynamics of fluorescent dyes and neuropeptides in the brain extracellular space. The models are encoded in Python using the reaction-diffusion module of the NEURON simulator and its recently developed extracellular reaction-diffusion simulator ([McDougal et al. 2013](https://doi.org/10.3389/fninf.2013.00028), [Newton et al. 2018](https://doi.org/10.3389/fninf.2018.00041)).

----

### Table of Contents

  1. [Install](#install)
       1. [Dependencies](#dependencies)
       2. [pip install](#pip-install)
  2. [The Models](#the-models)     
  3. [Usage details](#usage-details)
       1. [Importing](#importing)
       2. [Simulating](#simulating)
       3. [Chaning extracellular parameters](#changing-extracellular-parameters)       
  4. [License](#license)
  5. [Change Log](#change-log)
  6. [Citing](#citing)

----

## Install

**extracellular-models** installs as the `extmodels` Python package. It is tested with Python 3.9.

### Dependencies
Note that `extmodels` has the following core dependencies:
   * [NEURON](https://www.neuron.yale.edu/neuron)
   * [NumPy](http://www.numpy.org/)
   * [SciPy](https://www.scipy.org/)
   * [Matplotlib](https://matplotlib.org/)

### pip install
You can install `extmodels` with `pip` sourced from the GitHub repo:
```
pip install -e git+https://github.com/NTBEL/extracellular-models#egg=extmodels
```

----

## The Models
extracellular-models currently contains 3 models:
  * `extmodels.calcein` - Model of calcein diffusion in the brain extracellular space after photorelease from gold-coated nanovesicles.
  * `extmodels.somatastatin` - Model of somatastatin diffusion and loss (degradation, clearance, etc.) in the brain extracellular space after photorelease from gold-coated nanovesicles.
  * `extmodels.dynophin` - Model of dynorphin A diffusion and binding to the kLight receptor-based fluorescent sensor in the brain extracellular space after photorelease from gold-coated nanovesicles.     

----

## Usage details

To facilitate reusability and allow models to be importable they have been encapsulated as Python class objects and defined within separate modules.

### Importing
Models can be imported from their respective modules. For example,
```python
from extmodels.calcein import model
```
will import the initialized instance of the calcein model.

### Simulating
 To simulate the model you can call the `simulate` function. For example,
```python
time_step = 10 # time step in milliseconds
n_steps = 100 # run for 100 time steps.
out_frequency = 10 # compute model observables every 10 time steps.
model.simulate(time_step, n_steps, out_frequency)
```
will run the model simulation for total simulated time of 1 second (`time_step*n_step`) and store model observables at intervals of every 10 time steps. The observables and corresponding simulated times are accessible through the properteries `observables` and `times`; for example:
```python
model_observables = model.observables
times = model.times
```

### Changing Extracellular parameters
The NEURON rxd.Extracellular object that is used to define the extracellular space and its parameters (e.g., volume fraction and tortuosity) currently doesn't support modifying the parameters (such as volume fraction or tortuosity) after it has been initialized in the model. If you want to use a particular model but want, for example, a different tortuosity, you can import the `Model` class and then create a custom instance:
```python
from extmodels.calcein import Model

# Initialize the calcein model but with a tortuosity of 2 instead of the default
# 1.7 used for extmodels.calcein.model.
model = Model(tortuosity=2.0)
```

----

## Contact

Please open a [GitHub Issue](https://github.com/NTBEL/extracellular-models/issues) to
report any problems/bugs or make any comments, suggestions, or feature requests.

------

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for detail

----

## Change Log

See: [CHANGELOG](CHANGELOG)

----

## Citing
If any of these models are useful in your research, please cite this GitHub repo: https://github.com/NTBEL/extracellular-models
