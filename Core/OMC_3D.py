# OpenMC 0.13.3 Parameterised 3D Reactor Core with Reflector Model


############### Library Imports
from copy import deepcopy
from math import pi, sqrt
import openmc
import openmc.model


############## Define Model Builder Function
def OMC_model(MID : int = 0, RID : int = 3,
              MFR : float = 1, HLR : int = 2,
              RRT : float = 1, RHR : float = 0.5413) -> \
    tuple[openmc.Model, str, str, float, float, float, float]:
    ############### Material Definitions
    ######## Import Base Materials
    BaseMat = openmc.Materials.from_xml('BaseMat')

    # 19.75 % Enriched HALEU Uranium Nitride with 99.5 % N15 at 1600 K
    Fuel : openmc.Material = deepcopy(BaseMat[7])
    Fuel.depletable = True
    Fuel.temperature = 1600 # type: ignore

    # MA956 ODS Steel Clad Material at 1500 K
    Clad : openmc.Material = deepcopy(BaseMat[8])
    Clad.temperature = 1500 # type: ignore

    # 40 g/mol He Coolant at 1200 K and 2 MPa
    Gas : openmc.Material = deepcopy(BaseMat[9])
    Gas.temperature = 1200 # type: ignore
    Gas.set_density('g/cc', (1E-3 * 2E6) / ((8.314462 / 0.040) * Gas.temperature))

    # Moderator at 1200 K
    Mod : openmc.Material = deepcopy(BaseMat[MID])
    Mod.temperature = 1200 # type: ignore

    # Reflector at 900 K
    Ref : openmc.Material = deepcopy(BaseMat[RID])
    Ref.id = 99
    Ref.temperature = 900 # type: ignore


    ############### Geometry and Cell Definitions
    ######## Fuel Pin Geometry Definition
    # 1.58 cm Diameter Fuel Pellet
    fp_dia = 1.58 # cm
    s_fp = openmc.ZCylinder(0, 0, fp_dia / 2)
    c_fp = openmc.Cell(name='UN Fuel Pellet', fill=Fuel, region=-s_fp)

    # 1.68 cm Diameter Fuel Cladding
    cl_dia = 1.68 # cm
    s_cl = openmc.ZCylinder(0, 0, cl_dia / 2)
    c_cl = openmc.Cell(name='MA956 ODS Steel Cladding', fill=Clad, region=-s_cl & +s_fp)

    # 0.75 cm^2 HeXe Coolant Channel
    cc_dia = 2 * sqrt((0.75 / pi) + (cl_dia / 2)**2)
    s_cc = openmc.ZCylinder(0, 0, cc_dia / 2)
    c_cc = openmc.Cell(name='HeXe Coolant Channel', fill=Gas, region=-s_cc & +s_cl)

    # Moderator Region at Outside HeXe Coolant Channel
    c_mr = openmc.Cell(name='Moderator Region', fill=Mod, region=+s_cc)

    # Create Fuel Element Universe
    FE_univ = openmc.Universe(name='FE', cells=[c_fp, c_cl, c_cc, c_mr])

    ######## FE Hexagonal Lattice Definition
    # Create Hexagonal Lattice at Origin
    hex_lat = openmc.HexLattice(name='FE Lattice')
    hex_lat.center = (0.0, 0.0)

    # Create Cell with Reflector to fill Outskirts of the Lattice
    hex_lat.outer = openmc.Universe(cells=[openmc.Cell(name='Outer Reflector Sheath',
                                                       fill=Ref)])

    # Assign FE Universes to Lattice Positions from the Outermost Ring Inward
    # See Generating Function for Centered Hexagonal Numbers
    # https://en.wikipedia.org/wiki/Centered_hexagonal_number
    # H(n + 1) - H(n) = 6 * n
    # Use Input HLR which is the Number of Rings in the Hexagonal Lattice
    lat_univ = []
    i = HLR
    while i > 1:
        lat_univ.append([FE_univ] * (i - 1) * 6) # type: ignore
        i -= 1
    # Append Centermost FE Universe
    lat_univ.append([FE_univ])
    hex_lat.universes = lat_univ

    # Calculate and Assign the Lattice Pitch
    # Use Input MFR to Find Hexagonal Prism Edge Length
    a_hp = sqrt(((MFR * 0.25 * pi * fp_dia**2) + 0.25 * pi * cc_dia**2)
                       * (2 / (3 * sqrt(3))))

    # Use Edge Length to Determine the Circumscribing Circle's Diameter
    # This Value is Equivalent to the 2D Lattice Pitch
    hp_dia = 2 * a_hp
    hex_lat.pitch = (hp_dia,)

    # Calculate Radius of Circumscribing Circle of the Complete Core
    # Use Input RRT which is the Reflector Thickness
    # Use Input HLR which is the Number of Rings in the Hexagonal Lattice
    r_core : float = (HLR - 0.5) * hp_dia + RRT

    # Create Cylinder to Contain Hexagonal Lattice
    s_hex = openmc.ZCylinder(0, 0, r_core - RRT)

    # Create Radial Reflector Slab
    s_core = openmc.ZCylinder(0, 0, r_core, boundary_type='vacuum')

    # Define Z Planes to Limit Axial Extent of Core
    # Use Input RRT which is the Reflector Thickness
    # Use Input RHR which is the R/H Ratio
    # See Read 2020 for More Information on Default Value
    h_core : float = (r_core - RRT) / RHR
    s_top = openmc.ZPlane(h_core * 0.5, boundary_type='vacuum')
    s_bot = openmc.ZPlane(h_core * -0.5, boundary_type='vacuum')

    # Define 3D Cells for FE Lattice and Radial Reflector
    c_core = openmc.Cell(name='Core Lattice', fill=hex_lat,
                         region=-s_hex & -s_top & +s_bot)
    c_ref = openmc.Cell(name='Radial Reflector', fill=Ref,
                        region=-s_core & +s_hex & -s_top & +s_bot)

    ######## Define Root Universe of Geometry
    root_univ = openmc.Universe(name='Root Universe', cells=[c_core, c_ref])


    ############### Calculate Reactor Core Mass
    # Calculate Mass of Individual Fuel Element
    # Use Input MFR to Find Hexagonal Prism Area
    fe_mass = 0
    # Fuel Pellet Mass (g)
    fe_mass += 0.25 * pi * fp_dia**2 * h_core * Fuel.get_mass_density()
    # Cladding Mass (g) 
    fe_mass += 0.25 * pi * (cl_dia**2 - fp_dia**2) * h_core * Clad.get_mass_density()
    # Moderator Mass (g)
    fe_mass += (MFR * 0.25 * pi * fp_dia**2) * h_core * Mod.get_mass_density()

    # Calculate Number of Fuel Elements Inside Reactor Core
    # See Generating Function for Centered Hexagonal Numbers
    # https://en.wikipedia.org/wiki/Centered_hexagonal_number
    # Use Input HLR to Find the Number of Fuel Elements in the Hexagonal Lattice
    num_fe = ((3 * (HLR**2)) - (3 * HLR) + 1)
    hex_mass = num_fe * fe_mass

    # Calculate Mass of Radial Reflector (g)
    # Use Input MFR to Find Hexagonal Prism Area
    # Use Input RRT to Find Reflector Volume
    m_ref : float = (((pi * (r_core - RRT)**2\
                      - num_fe * (MFR * 0.25 * pi * fp_dia**2)\
                      + 0.25 * pi * cc_dia**2))\
                    + (pi * (2 * r_core * RRT - RRT**2))) * h_core\
                    * Ref.get_mass_density()
    
    # Calculate Total Core Mass in Kg
    m_core : float = (hex_mass + m_ref) * 1E-3

    # Convert Reflector Mass to Kg
    m_ref *= 1E-3


    ############### Initialise OpenMC Model
    model = openmc.Model()

    ######## Assign Materials and Geometry Data
    # Add Total Fuel Volume Data for Depletion
    Fuel.volume = num_fe * 0.25 * pi * fp_dia**2 * h_core # type: ignore
    model.materials = [Fuel, Clad, Gas, Mod, Ref]
    model.geometry = openmc.Geometry(root_univ)
    
    ######## Define Model Runtime Parameters
    # Setup Neutron Population and Criticality Cycle Parameters
    model.settings.batches = 30
    model.settings.inactive = 10
    model.settings.particles = 5000

    # Setup Cross Section Temperature Interpolation
    model.settings.temperature = {'method' : 'interpolation'}

    # Disable Tally Output Generation
    model.settings.output = {'tallies' : False}

    ######## Create Debug Geometry Plot
    # # XY Slice Plot
    # plt_dbgxy = openmc.Plot(name='Debug XY')
    # plt_dbgxy.basis = 'xy'
    # plt_dbgxy.color_by = 'material'
    # plt_dbgxy.colors = {Fuel : 'Yellow',
    #                     Clad : 'Grey',
    #                     Gas : 'Black',
    #                     Mod : 'Blue',
    #                     Ref : 'Red'}
    # plt_dbgxy.width = [100, 100]
    # plt_dbgxy.pixels = [8192, 8192]
    # plt_dbgxy.filename = 'XY_' + Mod.name + '_' + format(MFR, '.2f') \
    #                      + '_' + Ref.name + '_' + format(RRT, '.1f')

    # # XZ Slice Plot
    # plt_dbgxz = openmc.Plot(name='Debug XZ')
    # plt_dbgxz.basis = 'xz'
    # plt_dbgxz.color_by = 'material'
    # plt_dbgxz.colors = {Fuel : 'Yellow',
    #                     Clad : 'Grey',
    #                     Gas : 'Black',
    #                     Mod : 'Blue',
    #                     Ref : 'Red'}
    # plt_dbgxz.width = [100, 100]
    # plt_dbgxz.pixels = [8192, 8192]
    # plt_dbgxz.filename = 'XZ_' + Mod.name + '_' + format(MFR, '.2f') \
    #                      + '_' + Ref.name + '_' + format(RRT, '.1f')

    # # Plot the 3D Geometry
    # model.plots = [plt_dbgxy, plt_dbgxz]
    # model.plot_geometry()


    ############### Return OpenMC Model and Model Data
    return (model, Mod.name, Ref.name, r_core, h_core, m_ref, m_core)
