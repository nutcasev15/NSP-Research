# OpenMC 0.15.0 Parametrised Hexagonal Fuel Element Model


############### Library Imports
from math import pi, sqrt
import openmc
from openmc.model import HexagonalPrism, pin
from openmc.checkvalue import check_less_than, check_greater_than


############## Import Materials Database Symbols and Lookup Functions
from MatDB import *


############## Define Model Builder Function
def BuildModel(MID: int = 0, MFR : float = 1.0, TIT : float = 800.0) -> \
    tuple[openmc.Model, str, float, float, float, float]:

    ############### Reset OpenMC Object IDs for Materials
    openmc.reset_auto_ids()

    ############### Retrieve Material Definitions
    # 19.75 % Enriched HALEU Uranium Nitride with 99.5 % N15
    # Fuel Centreline Temperature Polynomial from Gas_Model.py
    FuelTemp : float = -7.6194E-13 * TIT**4 \
        + -1.1413E-08 * TIT**3 \
        + 7.6438E-05 * TIT**2 \
        + 8.1000E-01 * TIT \
        + 5.1174E+02
    Fuel : openmc.Material = get_mat_at_temp(MATID_UN, FuelTemp)
    Fuel.depletable = True

    # MA956 ODS Steel Clad Material
    # Clad Surface Temperature Polynomial from Gas_Model.py
    CladTemp : float = -8.1457E-12 * TIT**4 \
        + 3.2676E-08 * TIT**3 \
        + -3.4254E-05 * TIT**2 \
        + 9.6795E-01 * TIT \
        + 3.2131E+02
    Clad : openmc.Material = get_mat_at_temp(MATID_MA956ODS, CladTemp)

    # 40 g / mol HeXe Coolant at Turbine Inlet Temperature
    Gas : openmc.Material = get_mat_at_temp(MATID_HeXe, TIT)

    # Moderator at Turbine Inlet Temperature
    Mod : openmc.Material = get_mat_at_temp(MID, TIT)

    # Graphite Structural Material at Turbine Inlet Temperature
    Struct : openmc.Material = get_mat_at_temp(MATID_Graphite, TIT)

    ######## Define Moderator to Fuel Density Ratio
    MF_dens = Mod.get_mass_density() / Fuel.get_mass_density()


    ############### Geometry and Cell Definitions
    ######## Define Pellet Dimensions and Volumes
    # UN Fuel Pellet Diameter
    fp_dia = 1.58 # cm

    # Check Entered UN Fuel Pellet is Realistic with Tolerance
    # This Check is in SI Units
    check_less_than('UN Fuel Pellet Maximum Diameter',
                    fp_dia * 1E-2, 1.580E-2 + 0.01E-2)
    check_greater_than('UN Fuel Pellet Minimum Diameter',
                       fp_dia * 1E-2, 0.788E-2 - 0.01E-2)

    # Add Total Fuel Area Data to Fuel Material
    # Assign Area to Volume Field as this is a 2D Model
    Fuel.volume = 0.25 * pi * fp_dia**2

    # Calculate Moderator Pin Diameter for 6 Pin Arrangement
    # Use Input MFR to Find Moderator Pellet Area
    md_dia = sqrt(MFR * fp_dia**2 / 6) # cm

    # Check Calculated Moderator Pellet is Realistic with Tolerance
    # This Check is in SI Units
    check_less_than('Moderator Pellet Maximum Diameter',
                    md_dia * 1E-2, (1.580E-2 * sqrt(2.0)) + 0.01E-2)
    check_greater_than('Moderator Pellet Minimum Diameter',
                       md_dia * 1E-2, (0.788E-2 / sqrt(2.0)) - 0.01E-2)

    # Calculate Cladding Diameters
    cl_thick = 0.1 # cm
    fp_cl_dia = fp_dia + 2 * cl_thick # cm
    md_cl_dia = md_dia + 2 * cl_thick # cm

    # Calculate HeXe Fuel Coolant Channel Diameter
    fp_cc_area = 0.75 # cm^2
    fp_cc_dia = 2 * sqrt((fp_cc_area / pi) + (fp_cl_dia / 2)**2) # cm

    # Calculate Pellet Structural Sleeve Diameter
    slv_thick = 0.2 # cm
    fp_slv_dia = fp_cc_dia + 2 * slv_thick
    md_slv_dia = md_cl_dia + 2 * slv_thick


    ######## Define Pellet Surfaces
    s_fp = openmc.ZCylinder(0, 0, fp_dia / 2)
    s_fp_cl = openmc.ZCylinder(0, 0, fp_cl_dia / 2)
    s_fp_cc = openmc.ZCylinder(0, 0, fp_cc_dia / 2)
    s_fp_slv = openmc.ZCylinder(0, 0, fp_slv_dia / 2)

    s_md = openmc.ZCylinder(0, 0, md_dia / 2)
    s_md_cl = openmc.ZCylinder(0, 0, md_cl_dia / 2)
    s_md_slv = openmc.ZCylinder(0, 0, md_slv_dia / 2)


    ######## Create Fuel Pellet Universe
    u_fp = pin(surfaces=[s_fp, s_fp_cl, s_fp_cc, s_fp_slv],
               items=[Fuel, Clad, Gas, Struct, Struct])

    # Label Central Fuel Region in Fuel Pellet Universe
    cell_iter = iter(u_fp.cells)
    first_cell_id = next(cell_iter)
    u_fp.cells[first_cell_id].name = Fuel.name

    # Label Coolant Gas Region in Fuel Pellet Universe
    _ = next(cell_iter)
    third_cell_id = next(cell_iter)
    u_fp.cells[third_cell_id].name = Gas.name


    ######## Create Moderator Pellet Universe
    u_md = pin(surfaces=[s_md, s_md_cl, s_md_slv],
               items=[Mod, Clad, Struct, Struct])

    # Label Central Moderator Region in Moderator Pellet Universe
    first_cell_id = next(iter(u_md.cells))
    u_md.cells[first_cell_id].name = Mod.name


    ######## Moderator Pellet Hexagonal Lattice Definition
    # Create Hexagonal Lattice at Origin
    mod_lat = openmc.HexLattice(name='Moderator Pellet Lattice')
    mod_lat.orientation = 'x'
    mod_lat.center = (0.0, 0.0)

    # Create Cell with Structural Material to fill Lattice Outskirts
    mod_lat.outer = openmc.Universe(cells=[
        openmc.Cell(name='Outer Structural Sheath', fill=Struct)])

    # Assign Outer Lattice Positions to Moderator Pellet Universe
    lat_uni = []
    lat_uni.append([u_md] * 6)

    # Append Centremost Fuel Pellet Universe
    lat_uni.append([u_fp])
    mod_lat.universes = lat_uni

    # Use Edge Length to Determine the Circumscribing Circle's Diameter
    # This Value is Equivalent to the Moderator Lattice Pitch
    mod_lat.pitch = (max(fp_cc_dia, ((fp_slv_dia + md_slv_dia) / 2)), )

    ####### Hexagonal Carbonaceous Structural Boundary with Reflective BC
    # Calculate Circumscribing Radius of 2 Layer Moderator Pellet Lattice
    # Thickness Factor Ensures that the Edges are Tangent to the Pellets
    # Find Hexagonal Fuel Element Width Across the Flats
    lat_rad = (fp_slv_dia / 2) + md_slv_dia + slv_thick * 2
    hex_width = lat_rad * sqrt(3.0)
    s_cont = HexagonalPrism(lat_rad, 'x', (0, 0), 'reflective')

    # Fill Container with Moderator Pellet Lattice
    c_cont = openmc.Cell(name='Hexagonal Structural Clad', fill=mod_lat,
                         region=-s_cont)

    ####### Create Root Universe
    root_uni = openmc.Universe(name='Root Universe', cells=[c_cont])


    ############### Calculate Fuel Element Linear Mass
    fe_mass = 0.0

    ####### Calculate Fuel Pellet Pin Mass (g / cm)
    fe_mass += 0.25 * pi * fp_dia**2 * Fuel.get_mass_density()
    fe_mass += 0.25 * pi * (fp_cl_dia**2 - fp_dia**2) \
                    * Clad.get_mass_density()
    fe_mass += 0.25 * pi * (fp_slv_dia**2 - fp_cc_dia**2) \
                    * Struct.get_mass_density()

    ####### Calculate Moderator Pellet Pin Masses (g / cm)
    fe_mass += 6 * (0.25 * pi * md_dia**2 \
                        * Mod.get_mass_density() \
                    + 0.25 * pi * (md_cl_dia**2 - md_dia**2) \
                        * Clad.get_mass_density() \
                    + 0.25 * pi * (md_slv_dia**2 - md_cl_dia**2) \
                        * Struct.get_mass_density())

    # Add Hexagonal Structural Clad Mass (g / cm)
    # See https://en.wikipedia.org/wiki/Hexagon
    fe_mass += (0.8270 * pi * lat_rad**2 \
                - 0.25 * pi * fp_slv_dia**2 \
                - 6 * 0.25 * pi * md_slv_dia**2) \
                    * Struct.get_mass_density()

    # Convert Accumulated Mass to kg / m
    fe_mass *= 1E-1


    ############### Initialise OpenMC Model
    model = openmc.Model()

    ######## Assign Geometry Data
    model.geometry = openmc.Geometry(root_uni)

    ######## Assign Materials Data
    # Extract Subdivided Materials Data from Fuel Pellet Universe
    model.materials = list(u_fp.get_all_materials().values())
    # Add Materials Data for Moderator
    model.materials.extend([Mod])

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
    # plt_dbg.colors = {Fuel   : 'Green',
    #                   Clad   : 'Grey',
    #                   Struct : 'Black',
    #                   Gas    : 'Yellow',
    #                   Mod    : 'Blue'}
    # plt_dbg.width = [15, 15]
    # plt_dbg.pixels = [2440, 2440]
    # plt_dbg.filename = Mod.name + '_' + format(MFR, '.2f')
    # model.plots = [plt_dbg]
    # model.plot_geometry()


    ############### Return OpenMC Model and Model Data
    return (model, Mod.name, MF_dens, fp_cc_dia, hex_width, fe_mass, Mod, Fuel, Gas)
