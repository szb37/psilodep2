"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}

This repo is associated with this following publication
B. Szigeti et al.:
Assessing expectancy and suggestibility in a trial of escitalopram vs. psilocybin for depression treatment
Link: TBD

This script reproduces all figures of the paper.
For statistical models, see src/psilodep2 R project/.
"""

import src.folders as folders
import src.core as core
import pandas as pd
import os


df = pd.read_csv(os.path.join(folders.data, 'scores.csv'))

# Boxplots of pre-treatment variables - figure 1
if True:
    core.Plots.boxplots_per_treatment(
        df=df,
        variable='sss',
        y_label='Suggestibility')
    core.Plots.boxplots_per_treatment(
        df=df,
        variable='modtas',
        y_label='Absorption')
    core.Plots.boxplots_per_treatment(
        df=df,
        variable='rtx',
        y_label='Expectancy')

# Regression lines - figure 2 and 3
if True:
    core.Plots.regression(
        df=df,
        variable='rtx',
        xlim=[0, 100],
        x_label='Expectancy')
    core.Plots.regression(
        df=df,
        variable='sss',
        xlim=[25, 70],
        x_label='Suggestibility')
    core.Plots.regression(
        df=df,
        variable='modtas',
        xlim=[0, 80],
        x_label='Absorption')
