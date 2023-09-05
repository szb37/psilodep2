"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2023, B.Sz.
"""

import matplotlib.pyplot as plt
import seaborn as sns


scales = ['HAMD', 'BDI', 'MADRS', 'QIDS', 'STAIT', 'WEMWBS']

plt.rcParams.update({'font.family': 'arial'})
sns.set_style("darkgrid")
title_fontdict = {'fontsize': 18, 'fontweight': 'bold'}
axislabel_fontdict = {'fontsize': 16}
legend_fontsize = 11 #14
ticklabel_fontdict = {'fontsize': 14}
tick_labelsize = 14
