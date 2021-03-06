"""Deines the Base class for NEURON extracellular reaction-diffusion models.
"""

from abc import ABC, abstractmethod
from neuron import h, rxd
from neuron.units import nM, uM, cm, s, M, nm, um
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np


class ModelBase(ABC):
    """Abstract base class for extracellular reaction-diffusion models.

    Subclass this class when constructing new extracellular reaction-diffusion
    models.

    You can document each models components in its class docstring using the
    below example as a template.

    Model components (accessible as model attributes):

        Compartments:
            ecs -> rxd.Extracellular: The extracellular space domain and its
                parameters like volume fraction and tortuosity.
        Species:
            calcium -> rxd.Species: Ca2+ cation that diffuses in and leaks from
                the brain extracellular space.
        Parameters:
            loss_rate -> rxd.Parameter: Defines the kinetic rate for first order
                loss of calcium from the extracellular space.
        Reactions:
            ca_loss -> rxd.Rate: The reaction corresponding to first order loss
                of calcium from the extracellular space.
    """

    @abstractmethod
    def __init__(
        self,
        volume_fraction: float = 0.2,
        tortuosity: float = 1.7,
        dx: float = 5,
        xlo: float = -200.0,
        xhi: float = 200.0,
        ylo: float = -200.0,
        yhi: float = 200.0,
        zlo: float = -200.0,
        zhi: float = 200.0,
    ):
        """Use init to define and construct all the model components.

        You can follow the rxd Where/Who/What framework for model specification
        as in the example below. Make sure to also define private variables
        self._observables and self._times which will be used by
        the properties observables and times; here in init they can just be set
        to None but then they can used to store data in the simulate function.

        As in this example, you can make the key inputs for the rxd.Extracellular
        optional inputs to __init__ in case you want to create
        different versions of the model later with augmented ECS parameters, as
        NEURON currently doesn't support adjusting the rxd.Extracellular parameters
        (tortuosity, volume_fraction, etc.) after it has been initialized.

        Args:
            volume_fraction: The volume fraction of extracellular space
                (unitless). DEFAULT: 0.2
            tortuosity: The tortuosity for diffusion in the extracellular space
                (unitless). DEFAULT: 1.7
            dx: The discrettization size for volume voxels in the extracellular
                space in microns. DEFAUT: 5
            xlo: The position of the lower edge of the extracellular simulation
                domain along the x-direction in microns. DEFAULT: -200.
            xhi: The position of the upper edge of the extracellular simulation
                domain along the x-direction in microns. DEFAULT: 200.
            ylo: The position of the lower edge of the extracellular simulation
                domain along the y-direction in microns. DEFAULT: -200.
            yhi: The position of the upper edge of the extracellular simulation
                domain along the y-direction in microns. DEFAULT: -200.
            zlo: The position of the lower edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: -200.
            zhi: The position of the upper edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: -200.
        """
        # This is required for extracellular compartment simulations if using
        # NEURON version <= 7.6
        # rxd.options.enable.extracellular=True

        # Where? -- specify the regions
        # For extracellular reaction-diffusion we just have one
        # Extracellular region.
        self.ecs = rxd.Extracellular(
            xlo=xlo,
            ylo=ylo,
            zlo=zlo,
            xhi=xhi,
            yhi=yhi,
            zhi=zhi,
            dx=dx,
            volume_fraction=volume_fraction,
            tortuosity=tortuosity,
        )
        # Who? -- define all the species
        # Ca2+
        # The diffusion coefficient of Ca2+ is 5.3 x 10^-6 cm^2/s; https://doi.org/10.1016/0143-4160(87)90027-3
        d_calcium = 5.3e-6 * cm ** 2 / s  # diffusion coefficient
        ca_0 = 1 * uM  # initial concentration
        self.calcium = rxd.Species(
            self.ecs,
            name="calcium",
            d=d_calcium,
            charge=2,
            initial=ca_0,
            ecs_boundary_conditions=None,
            atolscale=1e-6,
        )
        # What? -- specify all the reactions
        # 1st-order loss of Ca2+
        kf = 1e-4 / s  # the reaction rate
        self.loss_rate = rxd.Parameter(ecs, value=kf)
        self.ca_loss = rxd.Rate(self.calcium, -self.loss_rate * self.calcium)

        # Initialize private variables -- required for all models
        self._times = None  # We'll store the time points in our simulation trajectory.
        self._observables = None  # We'll assign our observables to this variable.

    @abstractmethod
    def simulate(self, time_step: float, n_steps: int, output_frequency: int = 1):
        """Run the simulation and set the desired trajectory ouptuts.
        You can use the below example as a template.
        Args:
            time_step: The time interval between each step of the simulation
                trajectory in milliseconds. Times not in milliseconds can
                be converted using unit conversions from the neuron.units
                module.
            n_steps: The total number of simulation steps to run.
            output_frequency: Set the frequency for computing and storing
                any observables during the simulation. DEFAULT: 1

        """
        # Set the time step and initialize the simulator.
        h.dt = time_step
        h.finitialize()
        # define any observables that you want to store
        times = list()
        ca_zproject = list()
        # Run each step.
        for f in range(n_steps + 1):
            if (f == 0) or ((f % output_frequency) == 0):
                # simulation time
                times.append(f * time_step)
                # Get the mean value z-projection of the Ca2+ concentration.
                ca_zproject.append(self.calcium[self.ecs].states3d.mean(2))
            h.fadvance()
        self._times = np.array(times)
        self._observables = {"zproject_calcium": ca_zproject}
        return

    @staticmethod
    def animate_zprojection(
        zprojection_observable,
        min_value: float = 0,
        max_value: float = 1,
        interval: int = 10,
    ) -> animation.FuncAnimation:
        """Animates a z-projection observable over the simulation trajectory.
        Adapted from the runsim function in the Extracellular reaction-diffusion
        tutorial at:
            https://www.neuron.yale.edu/neuron/docs/extracellular-diffusion
        """
        frames = len(zprojection_observable)
        fig = plt.figure()
        im = plt.imshow(
            zprojection_observable[0], vmin=min_value, vmax=max_value, cmap="gray"
        )
        plt.axis("off")

        def init():
            im.set_data(zprojection_observable[0])
            return [im]

        def animate(i):
            im.set_data(zprojection_observable[i])
            return [im]

        anim = animation.FuncAnimation(
            fig, animate, init_func=init, frames=frames, interval=interval
        )
        return anim

    @property
    def observables(self):
        """Returns the dictionary of observables collected during simulation."""
        return self._observables

    @property
    def times(self):
        """Returns the time points in the simulation trajectory."""
        return self._times
