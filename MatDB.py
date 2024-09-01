# OpenMC 0.15.0 Model Base Material Database Script


############### Library Imports
import openmc
import openmc.checkvalue as cv
from numbers import Integral
import warnings


############### Define OpenMC Materials Container
MATDB = openmc.Materials()


############### Define Lookup Functions
def get_material_by_id(id) -> openmc.Material:
    cv.check_type('id', id, Integral)
    cv.check_greater_than('id', id, 0, equality=True)

    for mat in MATDB:
        if mat.id == id:
            return mat

    msg = f'No material instance exists with id={id}.'
    raise IndexError(msg)

def get_materials_by_name(name, case_sensitive=False, sort=False) \
    -> openmc.Materials:
    cv.check_type('name', name, str)

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


############### Base Material Definitions
# Graphite
Graph = openmc.Material(name='Graphite')
Graph.add_elements_from_formula('C')
Graph.set_density('g/cc', 1.710)
MATID_Graphite = Graph.id
MATDB.append(Graph)

# Zirconium Hydride
ZrH = openmc.Material(name='ZrH$_{1.66}$')
ZrH.add_elements_from_formula('Zr100H166')
ZrH.set_density('g/cc', 5.653)
MATID_ZrH = ZrH.id
MATDB.append(ZrH)

# Yttrium Hydride
YH = openmc.Material(name='YH$_{1.8}$')
YH.add_elements_from_formula('Y10H18')
YH.set_density('g/cc', 4.293)
MATID_YH = YH.id
MATDB.append(YH)

# Beryllium Oxide
BeO = openmc.Material(name='BeO')
BeO.add_elements_from_formula('BeO')
BeO.set_density('g/cc', 3.015)
MATID_BeO = BeO.id
MATDB.append(BeO)

# Magnesium Oxide
MgO = openmc.Material(name='MgO')
MgO.add_elements_from_formula('MgO')
MgO.set_density('g/cc', 3.581)
MATID_MgO = MgO.id
MATDB.append(MgO)

# Eutectic => BeO.MgO
BeO_MgO = openmc.Material.mix_materials([BeO, MgO], [0.7, 0.3], 'ao')
MATID_BeO_MgO = BeO_MgO.id
MATDB.append(BeO_MgO)

# Volumetric Composite => MgO-40ZrH
MgO_40ZrH = openmc.Material.mix_materials([MgO, ZrH], [0.6, 0.4], 'vo')
MgO_40ZrH.name = 'MgO-40ZrH$_{1.66}$'
MATID_MgO_40ZrH = MgO_40ZrH.id
MATDB.append(MgO_40ZrH)

# Volumetric Composite => MgO-40YH
MgO_40YH = openmc.Material.mix_materials([MgO, YH], [0.6, 0.4], 'vo')
MgO_40YH.name = 'MgO-40YH$_{1.8}$'
MATID_MgO_40YH = MgO_40YH.id
MATDB.append(MgO_40YH)

# Volumetric Composite => BeO.MgO-40ZrH
BeO_MgO_40ZrH = openmc.Material.mix_materials([BeO_MgO, ZrH], [0.6, 0.4], 'vo')
BeO_MgO_40ZrH.name = 'BeO.MgO-40ZrH$_{1.66}$'
MATID_BeO_MgO_40ZrH = BeO_MgO_40ZrH.id
MATDB.append(BeO_MgO_40ZrH)

# Volumetric Composite => BeO.MgO-40YH
BeO_MgO_40YH = openmc.Material.mix_materials([BeO_MgO, YH], [0.6, 0.4], 'vo')
BeO_MgO_40YH.name = 'BeO.MgO-40YH$_{1.8}$'
MATID_BeO_MgO_40YH = BeO_MgO_40YH.id
MATDB.append(BeO_MgO_40YH)

# Calcium Oxide
CaO = openmc.Material(name='CaO')
CaO.add_elements_from_formula('CaO')
CaO.set_density('g/cc', 3.340)
MATID_CaO = CaO.id
MATDB.append(CaO)

# Calcium Hydride
CaH = openmc.Material(name='CaH$_2$')
CaH.add_elements_from_formula('CaH2')
CaH.set_density('g/cc', 1.700)
MATID_CaH = CaH.id
MATDB.append(CaH)

# Volumetric Composite => CaO-40CaH
CaO_CaH = openmc.Material.mix_materials([CaO, CaH], [0.6, 0.4], 'vo')
CaO_CaH.name = 'CaO-40CaH$_2$'
MATID_CaO_CaH = CaO_CaH.id
MATDB.append(CaO_CaH)

# MA 956 ODS Steel
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

# 40 g/mol HeXe Gas
HeXe = openmc.Material(name='HeXe40')
HeXe.add_elements_from_formula('He71654Xe28346', 'ao')
# Gas has Temperature Dependant Density
MATID_HeXe = HeXe.id
MATDB.append(HeXe)

# 19.75 % Enriched HALEU Uranium Nitride with 99.5 % 15N
UN = openmc.Material(name='U$^{15}$N')
UN.add_nuclide('U235', 0.5 * 0.1975, 'ao')
UN.add_nuclide('U238', 0.5 * 0.8025, 'ao')
UN.add_element('N', 0.5, 'ao',
               enrichment=99.5,
               enrichment_target='N15',
               enrichment_type='ao')
UN.set_density('g/cc', 14.330)
MATID_UN = UN.id
MATDB.append(UN)
