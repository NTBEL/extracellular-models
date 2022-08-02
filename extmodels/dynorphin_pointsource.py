"""Model of dynorphin brain extracellular release, diffusion, and sensor-based detection.
"""

from neuron import h, rxd
from neuron.units import nM, uM, cm, s, M, nm, um
import numpy as np

# Avogadro's Number from scipy
from scipy.constants import N_A
from .basemodel import ModelBase


class Model(ModelBase):
    """A model of dynorphin A release, diffusion, and fluorescent sensor-based detection in the brain extracellular space.

    This model simulates the photorelease of the neuropeptide dynorphin A from
    gold coated nanovesicles, its subsequent diffusion in the brain
    extracellular space, and its interaction with the fluorescent receptor-based
    sensor kLight. An initial cylindrical bolus of dynorphin A is placed in
    the center of the extracellular domain to approximate the condition after
    photorelease from gold coated nanovesicles with release stimulated by a
    tornado scan of radius 30 micron with a multi-photon microscope. Dynorphin A
    can then diffuse away from the release site and can reversibly bind to the
    sensor protein which is assumed to have a uniform concentration in the
    brain extracellular space.

    Model components (accessible as model attributes):

        Compartments:
            ecs -> rxd.Extracellular: The extracellular space domain and its
                parameters like volume fraction and tortuosity.
        Species:
            dynorphin -> rxd.Species: Neuropeptide dynorphin A that is
                photoreleased in the brain extracellular space.
            sensor -> rxd.Species: The receptor-based sensor for dynophin.
            sensor_complex -> rxd.Species: The complex between dynorpin and the
                    sensor.
        Parameters:
            Kd -> rxd.Parameter: The dissociation constant for dynorphin-senson
                interaction.
            kon -> rxd.Parameter: The forward kinetic rate constant (on rate) for
                 second order binding of dynorphin to the sensor protein.
            koff -> rxd.Parameter: The reverse kinetic rate constant (off rate) for
                the first order unbinding of dynorphin from the sensor protein.
        Reactions:
            sensor_binding -> rxd.Rate: The reaction corresponding to reversible
                binding of dynorphin to the sensor protein.
                Forward reaction: dynorphin + sensor --kon--> dynorphin:sensor
                Reverse reaction: dynorphin:sensor --koff--> dynorphin + sensor

    """

    def __init__(
        self,
        volume_fraction: float = 0.2,
        tortuosity: float = 1.0,
        dx: float = 5,
        xlo: float = -202.5,
        xhi: float = 202.5,
        ylo: float = -202.5,
        yhi: float = 202.5,
        zlo: float = -202.5,
        zhi: float = 202.5,
    ):
        """Defines and constructs all the model components.

        Note:
            The simulation box edges should be defined such that there
            is an odd number of voxels along a given direction so that the
            central-most voxel will be centered at the origin (0,0,0) to get
            the expected point source initialization of dynorphin.

        Args:
            volume_fraction: The volume fraction of extracellular space
                (unitless). DEFAULT: 0.2
            tortuosity: The tortuosity for diffusion in the extracellular space
                (unitless). DEFAULT: 1.
            dx: The discretization size for volume voxels in the extracellular
                space in microns. DEFAUT: 5
            xlo: The position of the lower edge of the extracellular simulation
                domain along the x-direction in microns. DEFAULT: -202.5.
            xhi: The position of the upper edge of the extracellular simulation
                domain along the x-direction in microns. DEFAULT: 202.5.
            ylo: The position of the lower edge of the extracellular simulation
                domain along the y-direction in microns. DEFAULT: -202.5.
            yhi: The position of the upper edge of the extracellular simulation
                domain along the y-direction in microns. DEFAULT: 202.5.
            zlo: The position of the lower edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: -202.5.
            zhi: The position of the upper edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: 202.5.
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
        # Dynorphin (i.e., Dynorphin A(1-8))
        # The effective diffusion coefficeint for dynorphin is unknown, so we
        # we will assume something close to other similar peptides.
        d_dyn = (
            10.0e-7 * cm ** 2 / s
        )  # effective diffusion coefficient (factoring in tortuosity)

        # Estimated release is 1.2x10^8 molecs for SST:
        # Xiong et al. 2021 bioRxiv https://doi.org/10.1101/2021.09.10.459853
        # We can use the same for dynorphin.
        Qdyn = 1.2e8  # Number of molecules released during photostimulation

        # Approximate point-source release by initializing dynorphin in
        # the central-most voxel. Note the domain needs to be adjusted so
        # that the central voxel is centered on zero.
        voxel_volume = dx ** 3  # volume of one voxel
        dyn_0 = ((Qdyn / N_A) / voxel_volume) / volume_fraction * 1e15 * M
        dx_h = dx / 2
        self._dyn_0 = dyn_0
        self.dynorphin = rxd.Species(
            self.ecs,
            name="dynorphin A(1-8)",
            d=d_dyn,
            charge=0,
            initial=lambda nd: dyn_0
            if (nd.x3d ** 2 + nd.y3d ** 2 + nd.z3d ** 2 < dx_h ** 2)
            else 0,
            ecs_boundary_conditions=0.0,
        )
        # kLight - GPCR-based fluroescent sensor for dynophin
        # We'll assume that the sensor distribution can be approximated as
        # a uniform concentration throught the the brain ECS and that sensor
        # proteins don't diffuse.
        sensor_0 = 1 * nM  # sensor concentration
        self.sensor = rxd.Species(
            self.ecs, name="kLight", d=0, atolscale=1e-6, initial=sensor_0
        )
        # dynorphin bound to the kLight sensor
        self.sensor_complex = rxd.Species(
            self.ecs, name="kLight-dynorphin", d=0, initial=0, atolscale=1e-6
        )
        # What? -- specify all the reactions
        # dynorphin binding to the kLight sensor.
        # Dynorphin-receptor reaction radius - assume 2.5 nm, which is close
        # to half of the vertical height of KOR (~5 nm) roughly estimated from
        # the 6b73 PDB structure: https://www.rcsb.org/structure/6b73
        rad_on = 2.5e-7  # cm
        # Assume the binding rate is diffusion controlled:
        # 4piD*r => mL/s/molec. * 1e-3 => L/s/molec. * N_A (molec./mol) => 1/s/M
        kon = (
            (4 * np.pi * (d_dyn * s / cm ** 2) * rad_on * 1e-3) * N_A / (M * s)
        )  # on rate
        self.kon = rxd.Parameter(self.ecs, name="kon", value=kon)
        # Since we don't know the exact dissociation constant for
        # dynophin binding to the kLight sensor we can use the value
        # estimated for dynorphin binding to KOR which has Kd ~ 200 nM according to:
        # O'Connor et al. 2015 PNAS https://doi.org/10.1073/pnas.1510117112
        Kd_dyn = 200 * nM  # dissociation constant for dynorphin and kLight
        self.Kd = rxd.Parameter(self.ecs, name="Kd", value=Kd_dyn)
        self.koff = rxd.Parameter(
            self.ecs, name="koff", value=self.Kd.value * self.kon.value
        )  # off rate
        self.sensor_binding = rxd.Reaction(
            self.dynorphin + self.sensor, self.sensor_complex, self.kon, self.koff
        )

        # Initialize private variables -- required for all models
        self._times = None  # We'll store the time points in our simulation trajectory.
        self._observables = None  # We'll assign our observables to this variable.

    def simulate(
        self,
        time_step: float,
        n_steps: int,
        output_frequency: int = 1,
        nthread: int = 1,
    ):
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
        # Set the number of threads for the rxd module to use when simulating.
        rxd.nthread(nthread)
        # Set the time step and initialize the simulator.
        h.dt = time_step
        h.finitialize()
        # define any observables that you want to store
        times = list()
        dyn_zproject_mean = list()
        sensor_complex_zproject_mean = list()
        # Run each step.
        for f in range(n_steps + 1):
            if (f == 0) or ((f % output_frequency) == 0):
                # simulation time
                times.append(f * time_step)
                # Get the mean value z-projection of the dynorphin concentration.
                dyn_zproject_mean.append(self.dynorphin[self.ecs].states3d.mean(2))
                # Get the mean value z-projection of the dynorphin-sensor complex concentration.
                sensor_complex_zproject_mean.append(
                    self.sensor_complex[self.ecs].states3d.mean(2)
                )
            h.fadvance()
        self._times = np.array(times)
        self._observables = {
            "zproject_mean_dyn": dyn_zproject_mean,
            "zproject_mean_sensor_complex": sensor_complex_zproject_mean,
        }
        return

    @property
    def sensor_t_on(self):
        """Characteristic time for dynorphin binding to the kLight sensor."""
        return (1 / self.kon.value) / self._dyn_0

    @property
    def sensor_t_off(self):
        """Characteristic time for dynorphin dissociating from the kLight sensor."""
        return 1 / self.koff.value


model = Model()
