# OpenMC 0.13.3 XML Library Generator for Model Base Materials


############### Library Imports
import openmc


############### Base Material Definitions
# Graphite
Grph = openmc.Material(name='Graphite')
Grph.add_elements_from_formula('C')
Grph.set_density('g/cc', 1.710)
Grph.id = 0

# Zirconium Hydride
ZrH = openmc.Material(name='ZrH$_{1.66}$')
ZrH.add_elements_from_formula('Zr100H166')
ZrH.set_density('g/cc', 5.653)
ZrH.id = 1

# Yttrium Hydride
YH = openmc.Material(name='YH$_{1.8}$')
YH.add_elements_from_formula('Y10H18')
YH.set_density('g/cc', 4.293)
YH.id = 2

# Beryllium Oxide
BeO = openmc.Material(name='BeO')
BeO.add_elements_from_formula('BeO')
BeO.set_density('g/cc', 3.015)
BeO.id = 3

# Magnesium Oxide
MgO = openmc.Material(material_id=50, name='MgO')
MgO.add_elements_from_formula('MgO')
MgO.set_density('g/cc', 3.581)

# Volumetric Composite => MgO-40YH
MgO_40YH = openmc.Material.mix_materials([MgO, YH], [0.6, 0.4], 'vo')
MgO_40YH.id = 4
MgO_40YH.name = 'MgO-40YH$_{1.8}$'

# Eutectic => BeO.MgO
BeO_MgO = openmc.Material.mix_materials([BeO, MgO], [0.7, 0.3], 'ao')
BeO_MgO.id = 60

# Volumetric Composite => BeO.MgO-40YH
BeO_MgO_YH = openmc.Material.mix_materials([BeO_MgO, YH], [0.6, 0.4], 'vo')
BeO_MgO_YH.id = 5
BeO_MgO_YH.name = 'BeO.MgO-40YH$_{1.8}$'

# Calcium Oxide
CaO = openmc.Material(material_id=70, name='CaO')
CaO.add_elements_from_formula('CaO')
CaO.set_density('g/cc', 3.340)

# Calcium Hydride
CaH = openmc.Material(material_id=80, name='CaH$_2$')
CaH.add_elements_from_formula('CaH2')
CaH.set_density('g/cc', 1.700)

# Volumetric Composite => CaO-40CaH
CaO_CaH = openmc.Material.mix_materials([CaO, CaH], [0.6, 0.4], 'vo')
CaO_CaH.id = 6
CaO_CaH.name = 'CaO-40CaH$_2$'

# 19.75 % Enriched HALEU Uranium Nitride with 99.5 % 15N
UN = openmc.Material(name='U$^{15}$N')
UN.add_nuclide('U235', 0.5 * 0.1975, 'ao')
UN.add_nuclide('U238', 0.5 * 0.8025, 'ao')
UN.add_element('N', 0.5, 'ao',
               enrichment=99.5,
               enrichment_target='N15',
               enrichment_type='ao')
UN.set_density('g/cc', 14.330)

# MA 956 ODS Steel
MA = openmc.Material(name='MA956ODS')
MA.add_element('Fe', 74.45, 'wo')
MA.add_element('Cr', 20, 'wo')
MA.add_element('Al', 4.5, 'wo')
MA.add_element('Ti', 0.5, 'wo')
MA.add_element('Y', 0.5 * (2 / 5), 'wo')
MA.add_element('O', 0.5 * (3 / 5), 'wo')
MA.add_element('C', 0.05, 'wo')
MA.set_density('g/cc', 7.250)

# 40 g/mol HeXe Gas
HeXe = openmc.Material(name='HeXe40')
HeXe.add_elements_from_formula('He71654Xe28346', 'ao')
# Gas has Temperature Dependant Density

####### Export Materials Data to XML
openmc.Materials([Grph, ZrH, YH, BeO, MgO_40YH, BeO_MgO_YH, CaO_CaH,
                  UN, MA, HeXe]).export_to_xml('BaseMat')
