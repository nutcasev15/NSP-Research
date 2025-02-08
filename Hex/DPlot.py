# OpenMC 0.13.3 Python Post Processing Script for MF Ratio Depletion Study


############### Library Imports
import pickle
import itertools
import numpy as np
import matplotlib.markers as mrk
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt
import matplotlib.ticker as tck


############### Load Results from Pickle File
Pikname = 'Dep-Ref'
with open(Pikname + '.pkl', 'rb') as f:
    ModData = pickle.load(f)
    dep_sch = pickle.load(f)
    Mname = pickle.load(f)
    MF_den_ratio = pickle.load(f)
    channel_dia = pickle.load(f)
    hex_dia = pickle.load(f)
    hex_mass = pickle.load(f)
    K_EVL = pickle.load(f)
    K_BOL = pickle.load(f)
    K_EOL = pickle.load(f)
    U_burn = pickle.load(f)
    Coef_BOL = pickle.load(f)
    Coef_XST = pickle.load(f)
    Coef_MOL = pickle.load(f)
    Coef_EOL = pickle.load(f)


############### Plot Results
####### Define Set of Markers for Plotting
mrkrs = itertools.cycle(mrk.MarkerStyle.filled_markers)

###### Define Set of Line Styles for Plotting
lines = itertools.cycle(['-', '--', ':', '-.'])

####### Kinf Depletion History and Reactivity Coefficients
for (ModName, Kinf, BOL, XST, MOL, EOL) in zip(Mname, K_EVL, Coef_BOL, Coef_XST, Coef_MOL, Coef_EOL):
    fig = plt.figure(layout='constrained')
    ax = fig.add_subplot(111)
    ax.set_title('$K_{inf}$ Depletion History for ' + ModName)
    ax.set_xlabel('Burnup (MWd/kg)')
    ax.set_ylabel('$K_{inf}$')
    ax.grid()

    # Plot BOL Kinf Data with Error Bars for Each Moderator
    print(ModName + '\n')
    print(BOL)
    print('\n')
    # ax.errorbar(dep_sch, y=[j[0] for j in Kinf], yerr=[j[1] for j in Kinf])

    # # Finalise and Save BOL Kinf vs MF_vol Plot
    # fig.savefig(Pikname + '_' + ModName + '_Kinf.pdf', format='pdf')

    # # Clear Figure and Axes
    # fig.clear()
    # ax.clear()


####### Kinf and Burnup Data Table for Various Moderators
fig = plt.figure(figsize=(6, 1.5))
ax = fig.add_subplot(111)
ax.axis('off')
ax.grid('off')

# Retrieve Optimum Moderator to Fuel Volume Ratios
MF_vol = ['{:1.2f}'.format(T[1]) for T in ModData]

# Retrieve Kinf Values at BOL and EOL for Each Moderator
Kinf_BOL = ['{:1.3f}'.format(K[0][0]) for K in K_EVL]
Kinf_EOL = ['{:1.3f}'.format(K[-1][0]) for K in K_EVL]

# Calculate Difference in Kinf between BOL and EOL
K_diff = ['{:1.3f}'.format((K[0][0] - K[-1][0])) for K in K_EVL]

# Format Burnup Values for Table
Burnup = ['{:1.2f}'.format(B) for B in U_burn]

rows = ('Optimum MF Ratio', 'BOL $K_{inf}$', 'EOL $K_{inf}$', '$\\Delta K$', 'Burnup (U at. %)')

# Remove BeO Results and Assemble Cell Data
cell_text = [MF_vol[1:], Kinf_BOL[1:], Kinf_EOL[1:], K_diff[1:], Burnup[1:]]

Table_obj = ax.table(cellText=cell_text,
                     rowLabels=rows,
                     colLabels=Mname[1:],
                     loc='center', cellLoc='center')

Table_obj.auto_set_font_size(False)
Table_obj.set_fontsize(5.0)
for (row, col), cell in Table_obj.get_celld().items():
    if (row == 0 or col == -1):
        cell.set_text_props(fontproperties=fnt.FontProperties(weight='bold', size=3.5))

# Finalise and Save U Atom % Burnup vs MF_vol
# fig.legend(loc='outside right upper')
fig.savefig(Pikname + '_Table.pdf', format='pdf', bbox_inches='tight')

# Clear Figure and Axes
fig.clear()
ax.clear()
