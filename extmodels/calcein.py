"""Defines a model of calcein extracellular release and diffusion.
"""

from neuron import h, rxd
from neuron.units import nM, uM, cm, s, M, nm, um
import numpy as np
from .basemodel import ModelBase


class Model(ModelBase):
    """A model of calcein release and diffusion in the brain extracellular space.

    This model simulates the photorelease of the florescent dye calcein from
    gold coated nanovesicles and its subsequent diffusion in the brain
    extracellular space. An initial cylindrical bolus of calcein is placed in
    the center of the extracellular domain to approximate the condition after
    photorelease from gold coated nanovesicles with release stimulated by a
    tornado scan of radius 30 micron with a multi-photon microscope. Calcein
    can then diffuse away from the release site.

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
        # Approximate photorelease from gold coated nanovesicles similar to
        # SST in Xiong et al. 2021 bioRxiv https://doi.org/10.1101/2021.09.10.459853
        # Use an initial bolus of caclein with a uniform concentration inside a disc that
        # approximates the two-photon tornado scan area used for nanovesicle
        # release.
        r_stim = 30 * um  # radius of the tornado scan area
        cal_0 = 100 * uM  # concentration of calcein in the scan area after release.
        focal_disc_z = 12  # z-height of the 2-photon focal plane
        self.calcein = rxd.Species(
            self.ecs,
            name="calcein",
            d=d_calcein,
            charge=0,
            initial=lambda nd: cal_0
            if (nd.x3d ** 2 + nd.y3d ** 2 < r_stim ** 2)
            and (np.abs(nd.z3d) < focal_disc_z / 2)
            else 0,
            ecs_boundary_conditions=0.0,
        )
        # What? -- specify all the reactions
        # No reactions for calcein.

        # Initialize private variables -- required for all models
        self._times = None  # We'll store the time points in our simulation trajectory.
        self._observables = None  # We'll assign our observables to this variable.

    def simulate(self, time_step: float, n_steps: int, output_frequency: int = 1):
        """Run the simulation and set the desired trajectory ouptuts.

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
