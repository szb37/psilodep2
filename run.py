"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2020, DrugNerdsLab
:License: MIT

B. Szigeti et al.:
Assessing expectancy and suggestibility in a trial of escitalopram vs. psilocybin for depression treatment
Link: TBD

run.py:
  - fig 1, 2 and 3
  - supp figure 1

baseline_models.r:
  - supp table 1

models.r:
  - supp table 2-5; these tables contain all stats referenced in the paper
"""

import src.folders as folders
import src.core as core
import src.eqexp_analysis as eqexp
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

# Counterfactual analysis - what would be the results if expectancy would be equal, see Supplementary
# Counterfsctual data is generated using \psilodep2\codebase\src\R\get_equal_exp_data.r
if True:

    # Aggregate predicted outcomes data
    data_df = eqexp.EqExpCore.aggregate_predictions()

    # Calculate and aggregate model results across simulations
    stats_df = eqexp.EqExpCore.aggregate_model_results(
        data_df=pd.read_csv(os.path.join(folders.eqexp, 'predictions.csv')))

    # Get average results across iterations
    avg_df = eqexp.EqExpCore.get_avg_across_iterns(
        df=pd.read_csv(os.path.join(folders.eqexp, 'model_results.csv')))

    # Plot results
    eqexp.EqExpFigures.get_all_eqexp_figure(
        stats_df=pd.read_csv(os.path.join(folders.eqexp, 'model_results.csv')),
        avg_df=pd.read_csv(os.path.join(folders.eqexp, 'avg__model_results.csv')),)
