# OpenMC 0.15.2 Python Post Processing Script for MF Ratio Depletion Study


############### Library Imports
import pickle
import itertools
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt
import pandas as pd


############### Load Results from Pickle File
pkl_name = 'Dep-Ref'
with open(pkl_name + '.pkl', 'rb') as f:
  mod_data = pickle.load(f)
  dep_sch = pickle.load(f)
  mod_name = pickle.load(f)
  MF_den_ratio = pickle.load(f)
  channel_dia = pickle.load(f)
  hex_dia = pickle.load(f)
  hex_mass = pickle.load(f)
  K_EVL = pickle.load(f)
  K_BOL = pickle.load(f)
  K_EOL = pickle.load(f)
  U_burn = pickle.load(f)
  coef_BOL = pickle.load(f)
  coef_EOL = pickle.load(f)


############### Post Process Results
####### Define Set of Markers for Plotting
marks = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### Plot Kinf Depletion History
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('$K_{inf}$ Depletion History at 250 W/cm')
ax.set_xlabel('Operation Time (days)')
ax.set_xscale('log')
ax.set_ylabel('$K_{inf}$')
ax.grid()

for (name, Kinf) in zip(mod_name, K_EVL):
  # Plot BOL Kinf Data with Error Bars for Each Moderator
  ax.errorbar(
    x=dep_sch,
    y=[j[0] for j in Kinf],
    yerr=[j[1] for j in Kinf],
    label=name,
    marker=next(marks),
    linestyle=next(lines)
  )

# Finalise and Save Kinf vs Operation Time Plot
ax.legend(loc='lower left')
fig.tight_layout()
fig.savefig(pkl_name + '_Kinf.pdf', format='pdf')

# Clear Figure and Axes
fig.clear()
ax.clear()


####### Print Reactivity Coefficients
for (Name, BOL, EOL) in zip(mod_name, coef_BOL, coef_EOL):
  print(Name + '\n')
  print('BOL Temperature Reactivity Coefficients:')
  print(BOL)
  print('EOL Temperature Reactivity Coefficients:')
  print(EOL)
  print('\n')

####### Tabulate Depletion Summary for All Moderators
fig = plt.figure(figsize=(6, 1.5))
ax = fig.add_subplot(111)
ax.axis('off')
ax.grid('off')

# Retrieve Optimum Moderator to Fuel Volume Ratios
MF_vol = ['{:1.2f}'.format(T['MFR']) for T in mod_data]

# Format Optimum Linear Mass Value for Table
mass = ['{:1.3}'.format(M) for M in hex_mass]

# Retrieve Kinf Values at BOL and EOL for Each Moderator
Kinf_BOL = ['{:1.3f}'.format(K[0][0]) for K in K_EVL]
Kinf_EOL = ['{:1.3f}'.format(K[-1][0]) for K in K_EVL]

# Calculate Difference in Kinf between BOL and EOL
K_diff = ['{:1.3f}'.format((K[0][0] - K[-1][0])) for K in K_EVL]

# Format Burnup Values for Table
burnup = ['{:1.3f}'.format(B) for B in U_burn]

# Assign Row Labels for Table
rows = (
  'Optimum Volume Ratio',
  'Optimum Specific Mass (kg/m)',
  'BOL $K_{inf}$',
  'EOL $K_{inf}$',
  '$\\Delta K$',
  'Burnup (U at. $\\%$)'
)

# Assemble Cell Data
cell_text = [MF_vol, mass, Kinf_BOL, Kinf_EOL, K_diff, burnup]
Table_obj = ax.table(
  cellText=cell_text,
  rowLabels=rows,
  colLabels=mod_name,
  loc='center', cellLoc='center'
)

# Format Table Row and Column Labels
Table_obj.auto_set_font_size(False)
Table_obj.set_fontsize(5.0)
for (row, col), cell in Table_obj.get_celld().items():
  if (row == 0 or col == -1):
    cell.set_text_props(fontproperties=fnt.FontProperties(weight='bold', size=3.5))

# Finalise and Save Summary Table
fig.savefig(pkl_name + '_Table.pdf', format='pdf', bbox_inches='tight')

# Create Dataframe and Save Table as Plain LaTeX
pd.DataFrame(
  data=cell_text,
  index=rows,
  columns=mod_name
).to_latex(pkl_name + '_Table.tex')

# Clear Figure and Axes
fig.clear()
ax.clear()
