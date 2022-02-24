"""Defines a model of calcein extracellular release and diffusion.
"""

from neuron import h, rxd
from neuron.units import nM, uM, cm, s, M, nm, um
import numpy as np
# Avogadro's Number from scipy
from scipy.constants import N_A
from .basemodel import ModelBase


class Model(ModelBase):
    """A model of calcein release and diffusion in the brain extracellular space.

    This model simulates point source releas of the florescent dye calcein from
    and its subsequent diffusion in the brain extracellular space. An initial
    bolus of calcein the size of single simulation voxel is placed in the center
    of the extracellular domain to approximate an instantaneous point source.
    Calcein can then diffuse away from the source.

    Model components (accessible as model attributes):

        Compartments:
            ecs -> rxd.Extracellular: The extracellular space domain and its
                parameters like volume fraction and tortuosity.
        Species:
            calcein -> rxd.Species: A fluorescent dye that is photoreleased in
                the brain extracellular space.
        Parameters:
            None
        Reactions:
            None

    """

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
        """Defines and constructs all the model components.

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
                domain along the y-direction in microns. DEFAULT: 200.
            zlo: The position of the lower edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: -200.
            zhi: The position of the upper edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: 200.
        """
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
        # Calcein
        # Free diffusion from integrative optical imaging (IOI) is 44e-7 cm^2/s.
        d_calcein = 44e-7 * cm ** 2 / s  # diffusion coefficient
        # Approximate an instantaneous point source.
        Qcal = 1e8  # Number of calcein molecules released
        voxel_volume = dx**3 # volume of one voxel
        cal_0 = ((Qcal / N_A) / voxel_volume) / volume_fraction * 1e15 * M
        self.calcein = rxd.Species(
            self.ecs,
            name="calcein",
            d=d_calcein,
            charge=0,
            initial=lambda nd: cal_0
            if (nd.x3d ** 2 + nd.y3d ** 2 + nd.z3d**2 < dx ** 2) else 0,
            ecs_boundary_conditions=0.0,
        )
        # What? -- specify all the reactions
        # No reactions for calcein.

        # Initialize private variables -- required for all models
        self._times = None  # We'll store the time points in our simulation trajectory.
        self._observables = None  # We'll assign our observables to this variable.

    def simulate(self, time_step: float, n_steps: int,
                output_frequency: int = 1, nthread: int = 1,):
        """Run the simulation and set the desired trajectory ouptuts.

        Args:
            time_step: The time interval between each step of the simulation
                trajectory in milliseconds. Times not in milliseconds can
                be converted using unit conversions from the neuron.units
                module.
            n_steps: The total number of simulation steps to run.
            output_frequency: Set the frequency for computing and storing
                any observables during the simulation. DEFAULT: 1
            nthread: The number threads to use for multithreaded parallelization
                when running the simulation. DEFAULT: 1

        """
        # Set the number of threads for the rxd module to use when simulating.
        rxd.nthread(nthread)
        # Set the time step and initialize the simulator.
        h.dt = time_step
        h.finitialize()
        # define any observables that you want to store
        times = list()
        cal_zproject_mean = list()
        # Run each step.
        for f in range(n_steps + 1):
            if (f == 0) or ((f % output_frequency) == 0):
                # simulation time
                times.append(f * time_step)
                # Get the mean value z-projection of the Calcein concentration.
                cal_zproject_mean.append(self.calcein[self.ecs].states3d.mean(2))
            h.fadvance()
        self._times = np.array(times)
        self._observables = {"zproject_mean_calcium": cal_zproject_mean}
        return


model = Model()
