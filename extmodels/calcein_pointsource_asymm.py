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
        tortuosity_neg: float = 1.9,
        tortuosity_pos: float = 1.3,
        dx: float = 10,
        xlo: float = -405.0,
        xhi: float = 405.0,
        ylo: float = -405.0,
        yhi: float = 405.0,
        zlo: float = -405.0,
        zhi: float = 405.0,
        Qcal: int = 1e8,
    ):
        """Defines and constructs all the model components.

        Note: The simulation box edges should be defined such that there
        is an odd number of voxels along a given direction so that the
        central-most voxel will be centered at the origin (0,0,0) to get
        the expected point source initialization of calcien.

        Args:
            volume_fraction: The volume fraction of extracellular space
                (unitless). DEFAULT: 0.2
            tortuosity_neg: The tortuosity for diffusion in the extracellular space
                on the negative x half of the simulation box.
                (unitless). DEFAULT: 1.3
            tortuosity_pos: The tortuosity for diffusion in the extracellular space
                on the positive x half of the simulation box.
                (unitless). DEFAULT: 1.9
            dx: The discretization size for volume voxels in the extracellular
                space in microns. DEFAUT: 10
            xlo: The position of the lower edge of the extracellular simulation
                domain along the x-direction in microns. DEFAULT: -405.
            xhi: The position of the upper edge of the extracellular simulation
                domain along the x-direction in microns. DEFAULT: 405.
            ylo: The position of the lower edge of the extracellular simulation
                domain along the y-direction in microns. DEFAULT: -405.
            yhi: The position of the upper edge of the extracellular simulation
                domain along the y-direction in microns. DEFAULT: 405.
            zlo: The position of the lower edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: -405.
            zhi: The position of the upper edge of the extracellular simulation
                domain along the z-direction in microns. DEFAULT: 405.
            Qcal: The quantum of calcein molecules released from the point
                source.
        """

        # Define the voxel mesh here so we set the tortuosity value to be
        # different in different spatial regions of the simulation box.
        n_pixels = (np.array([xhi - xlo, yhi - ylo, zhi - zlo]) / dx).astype(int)
        xedges = np.linspace(xlo, xhi, n_pixels[0] + 1, endpoint=True)
        xpos = xedges[0:-1] + (xedges[1:] - xedges[:-1]) / 2
        yedges = np.linspace(ylo, yhi, n_pixels[1] + 1, endpoint=True)
        ypos = yedges[0:-1] + (yedges[1:] - yedges[:-1]) / 2
        zedges = np.linspace(zlo, zhi, n_pixels[2] + 1, endpoint=True)
        zpos = zedges[0:-1] + (zedges[1:] - zedges[:-1]) / 2
        # Generate a meshgrid for the x and y positions of pixels in micron
        xv, yv, zv = np.meshgrid(xpos, ypos, zpos, indexing="ij")
        self._xmesh = xv
        dx_h = dx / 2
        self._pos_xmask = xv > dx_h
        self._neg_xmask = xv < -dx_h
        tort = np.zeros_like(xv)
        tort[:, :, :] = tortuosity_neg
        tort[self._pos_xmask] = tortuosity_pos
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
            tortuosity=tort,
        )
        self._volume_fraction = volume_fraction
        # Who? -- define all the species
        # Calcein
        # Free diffusion from integrative optical imaging (IOI) is 44e-7 cm^2/s.
        d_calcein = 44e-7 * cm ** 2 / s  # diffusion coefficient
        # Approximate an instantaneous point source.
        voxel_volume = dx ** 3  # volume of one voxel
        cal_0 = ((Qcal / N_A) / voxel_volume) / volume_fraction * 1e15 * M
        self._voxel_volume = voxel_volume
        self.calcein = rxd.Species(
            self.ecs,
            name="calcein",
            d=d_calcein,
            charge=0,
            initial=lambda nd: cal_0
            if (nd.x3d ** 2 + nd.y3d ** 2 + nd.z3d ** 2 < dx_h ** 2)
            else 0,
            ecs_boundary_conditions=0.0,
        )
        # What? -- specify all the reactions
        # No reactions for calcein.

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
        cal_zproject_sum = list()
        # Number of calcein molecules that diffuse in negative x-direction.
        molec_neg = list()
        # Number of calcein molecules that diffuse in positive x-direction.
        molec_pos = list()
        microM_to_molec = (
            self._volume_fraction * 1e-6 * self._voxel_volume * 1e-15 * N_A
        )
        # Run each step.
        for f in range(n_steps + 1):
            if (f == 0) or ((f % output_frequency) == 0):
                # simulation time
                times.append(f * time_step)
                # Get the mean value z-projection of the Calcein concentration.
                cal_zproject_mean.append(self.calcein[self.ecs].states3d.mean(2))
                cal_zproject_sum.append(self.calcein[self.ecs].states3d.sum(2))
                num_neg = np.sum(
                    self.calcein[self.ecs].states3d[self._neg_xmask] * microM_to_molec
                )
                molec_neg.append(num_neg)
                num_pos = np.sum(
                    self.calcein[self.ecs].states3d[self._pos_xmask] * microM_to_molec
                )
                molec_pos.append(num_pos)
            h.fadvance()
        self._times = np.array(times)
        self._observables = {
            "zproject_mean_calcein": cal_zproject_mean,
            "zproject_sum_calcein": cal_zproject_sum,
            "num_neg": np.array(molec_neg),
            "num_pos": np.array(molec_pos),
        }
        return


model = Model()
