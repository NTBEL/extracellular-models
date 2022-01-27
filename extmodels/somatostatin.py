"""Model of somatostatin brain extracellular release, diffusion, and loss.
"""

from neuron import h, rxd
from neuron.units import nM, uM, cm, s, M, nm, um
import numpy as np

# Avogadro's Number from scipy
from scipy.constants import N_A
from .basemodel import ModelBase


class Model(ModelBase):
    """A model of somatostatin release, diffusion, and loss in the brain extracellular space.

    This model simulates the photorelease of the neuropeptide somatostatin from
    gold coated nanovesicles and its subsequent diffusion in the brain
    extracellular space. An initial cylindrical bolus of somatostatin is placed in
    the center of the extracellular domain to approximate the condition after
    photorelease from gold coated nanovesicles with release stimulated by a
    tornado scan of radius 30 micron with a multi-photon microscope. Somatostatin
    can then diffuse away from the release site and can be lost from the
    extracellular space through a first order decay reaction.

    Model components (accessible as model attributes):

        Compartments:
            ecs -> rxd.Extracellular: The extracellular space domain and its
                parameters like volume fraction and tortuosity.
        Species:
            sst -> rxd.Species: Neuropeptide somatostatin (sst) that is
                photoreleased in the brain extracellular space.
        Parameters:
            loss_rate -> rxd.Parameter: Defines the kinetic rate for first order
                loss of sst from the extracellular space.
        Reactions:
            sst_loss -> rxd.Rate: The reaction corresponding to first order loss
                of sst from the extracellular space.
                Forward reaction: sst --loss_rate--> None

    """

    def __init__(
        self,
        volume_fraction: float = 0.2,
        tortuosity: float = 1.0,
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
                (unitless). DEFAULT: 1.
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
        # The effective diffusion coefficient for fluorescently labelled SST
        # in acute brain slices from integrative optical imaging (IOI) is
        # 8.9+-1.3 x10^-7 cm^2/s. Xiong et al. 2021 bioRxiv https://doi.org/10.1101/2021.09.10.459853
        d_sst = (
            8.9e-7 * cm ** 2 / s
        )  # effective diffusion coefficient (factoring in tortuosity)
        # Approximate SST photorelease from gold coated nanovesicles as in
        # Xiong et al. 2021 bioRxiv https://doi.org/10.1101/2021.09.10.459853
        # Use an initial bolus of SST with a uniform concentration corresponding
        # to the estimated 1.2x10^8 released molecs inside a disc that
        # approximates the two-photon tornado scan area used for nanovesicle
        # release.
        r_stim = 30 * um  # radius of the tornado scan area
        Qsst = 1.2e8  # Number of SST molecules released during photostimulation
        focal_disc_z = 12  # z-height of the 2-photon focal plane
        disc_vol = (
            np.pi * r_stim ** 2 * focal_disc_z
        )  # volume of the tornado scan area.
        # SST concentration accounting for the disc volume and volume faction.
        sst_0 = ((Qsst / N_A) / disc_vol) / volume_fraction * 1e15 * M
        self.sst = rxd.Species(
            self.ecs,
            name="somatostatin",
            d=d_sst,
            charge=0,
            initial=lambda nd: sst_0
            if (nd.x3d ** 2 + nd.y3d ** 2 < r_stim ** 2)
            and (np.abs(nd.z3d) < focal_disc_z / 2)
            else 0,
            ecs_boundary_conditions=0.0,
        )
        # What? -- specify all the reactions
        # 1st-order loss of SST
        # Loss rate of SST estimated in range 0.023-0.048 per second.
        # Xiong et al. 2021 bioRxiv https://doi.org/10.1101/2021.09.10.459853
        kf = 3.6e-2 / s  # the reaction rate
        self.loss_rate = rxd.Parameter(self.ecs, value=kf)
        self.sst_loss = rxd.Rate(self.sst, -self.loss_rate * self.sst)

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
        sst_zproject_mean = list()
        # Run each step.
        for f in range(n_steps + 1):
            if (f == 0) or ((f % output_frequency) == 0):
                # simulation time
                times.append(f * time_step)
                # Get the mean value z-projection of the Calcein concentration.
                sst_zproject_mean.append(self.sst[self.ecs].states3d.mean(2))
            h.fadvance()
        self._times = np.array(times)
        self._observables = {"zproject_mean_sst": sst_zproject_mean}
        return


model = Model()
