"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import matplotlib.pyplot as pyplt
import seaborn as sns


scales = ['HAMD', 'BDI', 'MADRS', 'QIDS', 'STAIT', 'WEMWBS']
exp_types = ['avg', 'rtx']

pyplt.rcParams.update({'font.family': 'arial'})
sns.set_style("darkgrid")
title_fontdict = {'fontsize': 18, 'fontweight': 'bold'}
axislabel_fontdict = {'fontsize': 16}
legend_fontsize = 14
ticklabel_fontdict = {'fontsize': 14}
tick_labelsize = 14
