#TESTING SPM (This is for Paired TTest, will also test 2 sample)
import numpy as np                 
from matplotlib import pyplot 
import spm1d  
import pandas as pd
import matplotlib.pyplot as plt
from Hockey_Algo_Helper_File import *
from tkinter import Tk, filedialog
from collections import Counter
import matplotlib.cm as cm

SPM_DATA1 = filedialog.askopenfilename(title="Select Excel file", filetypes=[("Excel files", "*.xlsx;*.xls")])
VARA = pd.read_excel(SPM_DATA1, sheet_name='Z L5S1')
SPM_DATA2 = filedialog.askopenfilename(title="Select Excel file", filetypes=[("Excel files", "*.xlsx;*.xls")])
VARB = pd.read_excel(SPM_DATA2, sheet_name='Z L5S1')
#VARC = pd.read_excel(SPM_DATA, sheet_name='VariableC')

alpha      = 0.05
t          = spm1d.stats.ttest2(VARA, VARB, equal_var=True) #Parametric
ti         = t.inference(alpha, two_tailed=True, interp=True)
Joint = input("What side of body, what axis, and what joint is the Data from?: ")

#(2) Plot:
pyplot.close('all')
### plot mean and SD:
pyplot.figure( figsize=(8, 3.5) )
ax     = pyplot.axes( (0.1, 0.15, 0.35, 0.8) )
spm1d.plot.plot_mean_sd(VARA)
spm1d.plot.plot_mean_sd(VARB, linecolor='magenta', facecolor='pink')
ax.axhline(y=0, color='k', linestyle=':')
ax.set_xlabel('Shot Cycle (%)')
ax.set_ylabel(F'{Joint} (deg)', fontsize = 14)
# Adding colored bands
# Adding colored bands
#EDIT TO REFLECT THE EXACT SHOOTING SEGEMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.axvspan(0, 46, color='lightblue', alpha=1, ymin=0, ymax=0.05, label='BackSwing')
ax.axvspan(46, 75, color='darkblue', alpha=1, ymin=0, ymax=0.05, label= 'DownSwing' )
ax.axvspan(75, 100, color='red', alpha=1, ymin=0, ymax=0.05, label='Follow-Thru')


# Adding vertical lines
ax.axvline(x=46, color='lightblue', linestyle='--', label='Transition from Back Swing to Down Swing')
ax.axvline(x=75, color='red', linestyle='--', label='Point of Puck Contact')
ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')



############# plot SPM results: #########################################
ax     = pyplot.axes((0.55,0.15,0.35,0.8))
ti.plot()
ti.plot_threshold_label(fontsize=8)
ti.plot_p_values(size=10, offsets=[(0,0.3)])
ax.set_xlabel('Shot Cycle (%)')

ax.axvspan(0, 46, color='lightblue', alpha=1, ymin=0, ymax=0.05, label='BackSwing')
ax.axvspan(46, 75, color='darkblue', alpha=1, ymin=0, ymax=0.05, label= 'DownSwing' )
ax.axvspan(75, 100, color='red', alpha=1, ymin=0, ymax=0.05, label='Follow-Thru')

# Adding vertical lines
ax.axvline(x=46, color='lightblue', linestyle='--', label='Transition from Back Swing to Down Swing')
ax.axvline(x=75, color='red', linestyle='--', label='Point of Puck Contact')

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), shadow=False, ncol=2)  # Adjust the bbox_to_anchor and ncol as needed
# Create legend underneath the subplots
handles1, labels1 = ax.get_legend_handles_labels()
ax.legend(handles1, labels1, loc='upper center', bbox_to_anchor=(-0.25, -0.1), shadow=False, ncol=5)

pyplot.tight_layout()
pyplot.show()
pyplot.show()


