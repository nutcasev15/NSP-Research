# OpenMC 0.15.0 Parameterised Hexagonal Fuel Element Model


############### Library Imports
from copy import deepcopy
from math import pi, sqrt
import openmc
from openmc.model import HexagonalPrism


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############## Define Model Builder Function
def OMC_model(MID: int = 0, MFR : float = 1) -> \
    tuple[openmc.Model, str, float, float, float, float]:
    ############### Material Definitions
    # 19.75 % Enriched HALEU Uranium Nitride with 99.5 % N15 at 1600 K
    Fuel : openmc.Material = deepcopy(get_material_by_id(MATID_UN))
    Fuel.depletable = True
    Fuel.temperature = 1600 # type: ignore

    # MA956 ODS Steel Clad Material at 1500 K
    Clad : openmc.Material = deepcopy(get_material_by_id(MATID_MA956ODS))
    Clad.temperature = 1500 # type: ignore

    # 40 g/mol HeXe Coolant at 1200 K and 2 MPa
    Gas : openmc.Material = deepcopy(get_material_by_id(MATID_HeXe))
    Gas.temperature = 1200 # type: ignore
    Gas.set_density('g/cc', (1E-3 * 2E6) / ((8.314462 / 0.040) * Gas.temperature))

    # Moderator at 1200 K
    Mod : openmc.Material = deepcopy(get_material_by_id(MID))
    Mod.temperature = 1200 # type: ignore

    ######## Define Moderator to Fuel Density Ratio
    MF_dens : float = Mod.get_mass_density() / Fuel.get_mass_density()


    ############### Geometry and Cell Definitions
    # 1.58 cm Diameter Fuel Pellet with 3 Subdivisions for Depletion
    fp_dia = 1.58 # cm
    s_fp = openmc.ZCylinder(0, 0, fp_dia / 2)

    # Add Total Fuel Area Data to Fuel Material
    Fuel.volume = 0.25 * pi * fp_dia**2 # type: ignore

    # Create Fuel Pellet with 3 Subdivisions of Burnable Material Instances
    u_fp = openmc.model.pin([s_fp], [Fuel, Clad], {0 : 3}, divide_vols=True)
    c_fp = openmc.Cell(name='UN Fuel Pellet', fill=u_fp, region=-s_fp)

    # 1.68 cm Diameter Fuel Cladding
    cl_dia = 1.68 # cm
    s_cl = openmc.ZCylinder(0, 0, cl_dia / 2)
    c_cl = openmc.Cell(name='MA956 ODS Steel Cladding', fill=Clad, region=-s_cl & +s_fp)

    # 0.75 cm^2 HeXe Coolant Channel
    cc_dia : float = 2 * sqrt((0.75 / pi) + (cl_dia / 2)**2)
    s_cc = openmc.ZCylinder(0, 0, cc_dia / 2)
    c_cc = openmc.Cell(name='HeXe Coolant Channel', fill=Gas, region=-s_cc & +s_cl)

    # Hexagonal Prism Moderator at Origin with Reflective BC
    # Use Input MF Ratio to Find Hexagonal Prism Edge Length
    # Also Calculate Diameter of Circumscribing Circle
    a_hp = sqrt(((MFR * 0.25 * pi * fp_dia**2) + 0.25 * pi * cc_dia**2)
                * (2 / (3 * sqrt(3))))
    hp_dia : float = 2 * a_hp
    r_hp = HexagonalPrism(a_hp, 'y', (0, 0), 'reflective')
    c_mh = openmc.Cell(name='Hexagonal Moderator Block', fill=Mod, region=-r_hp & +s_cc)

    # Create Root Universe
    root_univ = openmc.Universe(name='Root Universe', cells=[c_fp, c_cl, c_cc, c_mh])


    ############### Calculate Fuel Element Linear Mass
    # Use Input MFR to Find Hexagonal Moderator Block Area
    fe_mass : float = 0
    # Fuel Pellet Mass (g/cm)
    fe_mass += 0.25 * pi * fp_dia**2 * Fuel.get_mass_density()
    # Cladding Mass (g/cm) 
    fe_mass += 0.25 * pi * (cl_dia**2 - fp_dia**2) * Clad.get_mass_density()
    # Moderator Mass (g/cm)
    fe_mass += (MFR * 0.25 * pi * fp_dia**2) * Mod.get_mass_density()
    # Convert to Kg/m
    fe_mass *= 1E-1


    ############### Initialise OpenMC Model
    model = openmc.Model()

    ######## Assign Geometry Data
    model.geometry = openmc.Geometry(root_univ)

    ######## Assign Materials Data
    # Extact Subdivided Materials Data from Fuel Pellet Universe
    model.materials = list(u_fp.get_all_materials().values())
    # Add Materials Data for Coolant and Moderator
    model.materials.extend([Gas, Mod])
    
    ######## Define Model Runtime Parameters
    # Setup Neutron Population and Criticality Cycle Parameters
    model.settings.batches = 30
    model.settings.inactive = 5
    model.settings.particles = 5000

    # Setup Cross Section Temperature Interpolation
    model.settings.temperature = {'method' : 'interpolation'}

    # Disable Tally Output Generation
    model.settings.output = {'tallies' : False}

    ######## Create Debug Geometry Plot
    # plt_dbg = openmc.Plot(name='Debug XY')
    # plt_dbg.basis = 'xy'
    # plt_dbg.color_by = 'material'
    # plt_dbg.colors = {Fuel : 'Yellow', Clad : 'Grey', Gas : 'Black', Mod : 'Blue'}
    # plt_dbg.width = [15, 15]
    # plt_dbg.pixels = [2440, 2440]
    # plt_dbg.filename = Mod.name + '_' + format(MFR, '.2f')
    # model.plots = [plt_dbg]
    # model.plot_geometry()


    ############### Return OpenMC Model and Model Data
    return (model, Mod.name, MF_dens, cc_dia, hp_dia, fe_mass)
