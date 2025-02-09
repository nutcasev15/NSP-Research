# OpenMC 0.15.0 Materials Database and Lookup Functions


############### Library Imports
import openmc
import openmc.checkvalue as cv
from numbers import Integral
import warnings


############### Define OpenMC Materials Container
MATDB = openmc.Materials()


############### Base Material Definitions
# Alpha Beryllium Oxide
# See https://doi.org/10.1016/j.net.2022.07.017
# Theoretical Density in Section 2.1
BeO = openmc.Material(name='BeO')
BeO.add_elements_from_formula('BeO')
BeO.set_density('g/cc', 3.010)
MATID_BeO = BeO.id
MATDB.append(BeO)

# Calcium Oxide
# See http://dx.doi.org/10.1155/2014/123478
# Table 3: CaO Density at 300 K
CaO = openmc.Material(name='CaO')
CaO.add_elements_from_formula('CaO')
CaO.set_density('g/cc', 3.350)
MATID_CaO = CaO.id
MATDB.append(CaO)

# Calcium Hydride [H/M = 2]
# See https://doi.org/10.1016/j.ijhydene.2023.04.088
# Table 1: Density at 298 K
CaH = openmc.Material(name='CaH$_2$')
CaH.add_elements_from_formula('CaH2')
CaH.set_density('g/cc', 1.700)
MATID_CaH = CaH.id
MATDB.append(CaH)

# 2D Woven Carbon Carbon Composite
# See http://dx.doi.org/10.1007/s10853-006-1123-3
# Table 1: 2D Woven Sample at 298 K
CCComp = openmc.Material(name='CarbonCarbon')
CCComp.add_elements_from_formula('C')
CCComp.set_density('g/cc', 1.720)
MATID_CarbonCarbon = CCComp.id
MATDB.append(CCComp)

# Nuclear Graphite GR-280
# See https://www-pub.iaea.org/MTCD/Publications/PDF/IAEA-THPH_web.pdf
# Table 4.3: Mean Density at 293 K
Graph = openmc.Material(name='Graphite')
Graph.add_elements_from_formula('C')
Graph.set_density('g/cc', 1.710)
MATID_Graphite = Graph.id
MATDB.append(Graph)

# 40 g/mol HeXe Gas
# Formula Calculated from Standard Atomic Weights and Molar Fractions
HeXe = openmc.Material(name='HeXe40')
HeXe.add_elements_from_formula('He71654Xe28346', 'ao')
MATID_HeXe = HeXe.id
MATDB.append(HeXe)

# MA 956 ODS Steel
# See https://doi.org/10.1016/j.jnucmat.2004.10.118
# Table 1: Elemental Composition
# Table 4: Density at 298 K
ODS = openmc.Material(name='MA956ODS')
ODS.add_element('Fe', 74.45, 'wo')
ODS.add_element('Cr', 20, 'wo')
ODS.add_element('Al', 4.5, 'wo')
ODS.add_element('Ti', 0.5, 'wo')
ODS.add_element('Y', 0.5 * (2 / 5), 'wo')
ODS.add_element('O', 0.5 * (3 / 5), 'wo')
ODS.add_element('C', 0.05, 'wo')
ODS.set_density('g/cc', 7.250)
MATID_MA956ODS = ODS.id
MATDB.append(ODS)

# Magnesium Oxide
# See http://dx.doi.org/10.1155/2014/123478
# Table 3: MgO Density at 300 K
MgO = openmc.Material(name='MgO')
MgO.add_elements_from_formula('MgO')
MgO.set_density('g/cc', 3.580)
MATID_MgO = MgO.id
MATDB.append(MgO)

# 19.75 % Enriched HALEU Uranium Nitride
# 99.5 % 15N Enriched Nitrogen for Cross Section Reduction
UN = openmc.Material(name='U$^{15}$N')
UN.add_nuclide('U235', 0.5 * 0.1975, 'ao')
UN.add_nuclide('U238', 0.5 * 0.8025, 'ao')
UN.add_element('N', 0.5, 'ao',
               enrichment=99.5,
               enrichment_target='N15',
               enrichment_type='ao')
MATID_UN = UN.id
MATDB.append(UN)

# Yttrium Hydride
# See https://www.osti.gov/biblio/2204161
# Table 8-2: Theoretical Density in Caption
# Table 8-3: Average H/M Ratio of Samples
YH = openmc.Material(name='YH$_{1.8}$')
YH.add_elements_from_formula('Y10H18')
YH.set_density('g/cc', 4.280)
MATID_YH = YH.id
MATDB.append(YH)

# Zirconium Hydride
# See https://doi.org/10.1016/S0925-8388(01)01448-7
# Table 1: Density at 273 K for H/M = 1.66
ZrH = openmc.Material(name='ZrH$_{1.66}$')
ZrH.add_elements_from_formula('Zr100H166')
ZrH.set_density('g/cc', 5.653)
MATID_ZrH = ZrH.id
MATDB.append(ZrH)


############### Mixed Material Definitions
# Eutectic BeO.MgO Matrix at 70 mol Percent BeO
# See https://doi.org/10.1111/j.1151-2916.1963.tb11708.x
BeO_MgO = openmc.Material.mix_materials([BeO, MgO], [0.7, 0.3], 'ao')
MATID_BeO_MgO = BeO_MgO.id
MATDB.append(BeO_MgO)

# Composite Neutron Moderator BeO.MgO-40YH
# See https://doi.org/10.1080/21870764.2021.1993592
# Maximum Entrained Moderator at 40 Volume Percent
BeO_MgO_40YH = openmc.Material.mix_materials([BeO_MgO, YH], [0.6, 0.4], 'vo')
BeO_MgO_40YH.name = 'BeO.MgO-40YH$_{1.8}$'
MATID_BeO_MgO_40YH = BeO_MgO_40YH.id
MATDB.append(BeO_MgO_40YH)

# Composite Neutron Moderator BeO.MgO-40ZrH
# See https://doi.org/10.1080/21870764.2021.1993592
# Maximum Entrained Moderator at 40 Volume Percent
BeO_MgO_40ZrH = openmc.Material.mix_materials([BeO_MgO, ZrH], [0.6, 0.4], 'vo')
BeO_MgO_40ZrH.name = 'BeO.MgO-40ZrH$_{1.66}$'
MATID_BeO_MgO_40ZrH = BeO_MgO_40ZrH.id
MATDB.append(BeO_MgO_40ZrH)

# Composite Neutron Moderator CaO-40CaH
# See https://doi.org/10.1007/s42452-020-03942-1 for Feasibility
# See https://doi.org/10.1016/0009-2509(89)85232-7 for Sintering
# Assumed 40 Volume Percent Entrainment is Possible
CaO_CaH = openmc.Material.mix_materials([CaO, CaH], [0.6, 0.4], 'vo')
CaO_CaH.name = 'CaO-40CaH$_2$'
MATID_CaO_CaH = CaO_CaH.id
MATDB.append(CaO_CaH)

# Composite Neutron Moderator MgO-40YH
# See https://doi.org/10.1080/21870764.2021.1993592
# Maximum Entrained Moderator at 40 Volume Percent
MgO_40YH = openmc.Material.mix_materials([MgO, YH], [0.6, 0.4], 'vo')
MgO_40YH.name = 'MgO-40YH$_{1.8}$'
MATID_MgO_40YH = MgO_40YH.id
MATDB.append(MgO_40YH)

# Composite Neutron Moderator MgO-40ZrH
# See https://doi.org/10.1080/21870764.2021.1993592
# Maximum Entrained Moderator at 40 Volume Percent
MgO_40ZrH = openmc.Material.mix_materials([MgO, ZrH], [0.6, 0.4], 'vo')
MgO_40ZrH.name = 'MgO-40ZrH$_{1.66}$'
MATID_MgO_40ZrH = MgO_40ZrH.id
MATDB.append(MgO_40ZrH)


############### Define Material Retrieval Function
def get_mat_at_temp(id : Integral, T : float) -> openmc.Material:
    ''' Retrieve Material from Database at Specified Temperature

        Clones the material instance in the database and
        modifies its density to account for the temperature.
        Also, updates the material nuclide temperature.

        Parameters
        ----------
        id : int
            The id to match in the database
        T : float
            The specified temperature of the material

        Returns
        -------
        openmc.Material
            Material of matching id with updated properties
            depending on the specified temperature

    '''
    cv.check_type('id', id, Integral)
    cv.check_greater_than('id', id, 0, equality=True)
    cv.check_greater_than('T', T, 0, equality=True)

    # Create New Instance of Material
    mat = get_material_by_id(id).clone()

    # Step Through All IDs to Update Temperature Dependent Properties
    # Worst Case Scenario Properties Assumed
    if id == MATID_BeO:
        # See https://doi.org/10.1016/j.net.2022.07.017
        # Figure 5: Max CTE is 14E-6 K^-1
        # Density Reference Temperature is 294 K
        rho_T = mat.get_mass_density() / (1 + 14E-6 * (T - 294.0))**3
    elif id == MATID_BeO_MgO or \
    id == MATID_BeO_MgO_40YH or \
    id == MATID_BeO_MgO_40ZrH:
        # Assuming Expansion Behaviour Similar to MgO
        # Assuming Matrix Expansion Behaviour is Dominant
        # See http://dx.doi.org/10.1155/2014/123478
        # Figure 6: Average CTE is 16E-6 K^-1 at 1600 K
        # Density Reference Temperature is 300 K
        rho_T = mat.get_mass_density() / (1 + 16E-6 * (T - 300.0))**3
    elif id == MATID_CaH:
        # No Density or CTE Data vs Temperature Available
        # Assuming No Changes
        rho_T = mat.get_mass_density()
    elif id == MATID_CaO or \
    id == MATID_CaO_CaH:
        # Assuming Matrix Expansion Behaviour is Dominant
        # See http://dx.doi.org/10.1155/2014/123478
        # Figure 6: Average CTE is 15E-6 K^-1 at 1600 K
        # Density Reference Temperature is 300 K
        rho_T = mat.get_mass_density() / (1 + 15E-6 * (T - 300.0))**3
    elif id == MATID_CarbonCarbon:
        # See http://dx.doi.org/10.1007/s10853-006-1123-3
        # Figure 3: Max CTE is 6E-6 K^-1
        # Density Reference Temperature is 298 K
        rho_T = mat.get_mass_density() / (1 + 6E-6 * (T - 298.0))**3
    elif id == MATID_Graphite:
        # See https://www-pub.iaea.org/MTCD/Publications/PDF/IAEA-THPH_web.pdf
        # Table 4.7: Max CTE is 6.3E-6
        # Density Reference Temperature is 293 K
        rho_T = mat.get_mass_density() / (1 + 6.3E-6 * (T - 293.0))**3
    elif id == MATID_HeXe:
        # Ideal Gas has Temperature Dependant Density
        # Density of 40 g/mol Ideal Gas at 2 MPa
        # Result Converted to g/cc
        rho_T = (1E-3 * (2E6 / ((8.314462 / 0.040) * T)))
    elif id == MATID_MA956ODS:
        # See https://doi.org/10.1016/j.jnucmat.2004.10.118
        # Table 4: Average CTE is 15.5E-6 at 1373 K
        # Density Reference Temperature is 298 K
        rho_T = mat.get_mass_density() / (1 + 15.5E-6 * (T - 298.0))**3
    elif id == MATID_MgO or \
    id == MATID_MgO_40YH or \
    id == MATID_MgO_40ZrH:
        # Assuming Matrix Expansion Behaviour is Dominant
        # See http://dx.doi.org/10.1155/2014/123478
        # Figure 6: Average CTE is 16E-6 K^-1 at 1600 K
        # Density Reference Temperature is 300 K
        rho_T = mat.get_mass_density() / (1 + 16E-6 * (T - 300.0))**3
    elif id == MATID_UN:
        # See https://www-pub.iaea.org/MTCD/Publications/PDF/IAEA-THPH_web.pdf
        # Polynomial Valid from 298 K to 2523 K
        # Result Converted to g/cc
        rho_T = (1E-3 * (14420 - 0.2779 * T - 4.897E-5 * T**2))
    elif id == MATID_YH:
        # See https://www.osti.gov/biblio/2204161
        # Figure 5-2: Max Instantaneous CTE is 25E-6
        # Density Reference Temperature is 300 K
        rho_T = mat.get_mass_density() / (1 + 25E-6 * (T - 300.0))**3
    elif id == MATID_ZrH:
        # See https://doi.org/10.1016/S0925-8388(01)01448-7
        # Table 1: Average CTE is 27.4E-6 at H/M = 1.66
        # Density Reference Temperature is 273 K
        rho_T = mat.get_mass_density() / (1 + 27.4E-6 * (T - 273.0))**3
    else:
        # Material not Found
        msg = f'No material exists in the database with id={id}.'
        raise IndexError(msg)

    # Assign Updated Material Properties
    mat.temperature = T
    mat.set_density('g/cc', rho_T)

    return mat


############### Define Database Lookup Functions
def get_material_by_id(id : Integral) -> openmc.Material:
    ''' Return Material from Database with Specified ID

        Parameters
        ----------
        id : int
            The id to match in the database

        Returns
        -------
        openmc.Material
            Material with matching id
    '''
    cv.check_type('id', id, Integral)
    cv.check_greater_than('id', id, 0, equality=True)

    for mat in MATDB:
        if mat.id == id:
            return mat

    msg = f'No material exists in the database with id={id}.'
    raise IndexError(msg)


def get_materials_by_name(name : str, case_sensitive : bool = False,
                          sort : bool = False) -> openmc.Materials:
    ''' Return a List of Materials from the Database with Matching Names

        Parameters
        ----------
        name : str
            The name to match
        case_sensitive : bool
            Whether to distinguish upper and lower case letters
            in each material's name. Letter case is not distinguished
            by default.
        sort : bool, optional
            Whether to sort the matched materials by id
            number to establish a predictable ordering

        Returns
        -------
        openmc.Materials
            A list of materials matching the queried name
    '''
    cv.check_type('name', name, str)
    cv.check_type('case_sensitive', case_sensitive, bool)
    cv.check_type('sort', sort, bool)

    matches = openmc.Materials()
    match_name = name if case_sensitive else name.lower()
    for mat in MATDB:
        mat_name = mat.name if case_sensitive else mat.name.lower()
        if mat_name == match_name:
            matches.append(mat)
    if sort is True:
        matches.sort(key=lambda x: x.id)

    if len(matches) == 0:
        msg = f'No materials found with names that match {match_name}.'
        warnings.warn(msg, UserWarning)

    return matches
