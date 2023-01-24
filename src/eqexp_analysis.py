"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import matplotlib.pyplot as pyplt
import src.folders as folders
from tqdm.contrib.itertools import product as tqdmproduct
import statsmodels.formula.api as smf
import src.config as config
from tqdm import tqdm
import itertools
import seaborn as sns
import pandas as pd
import numpy as np
import warnings
import copy
import os

# Some models likely will not converge in the counterfactual analysis, suppress warnings
warnings.filterwarnings("ignore")


class EqExpCore():

    def aggregate_predictions():
        ''' Generate single CSV containing all model results '''

        dfs = [pd.DataFrame(columns=['patient_id', 'trt', 'scale', 'tp',
                            'score', 'exp_value', 'exp_type', 'target_exp', 'itern'])]

        for filename in tqdm(os.listdir(folders.predictions), desc='Aggregate simulation data'):

            if 'predictions_' not in filename:
                continue

            scale, exp_type, target_exp, itern = EqExpCore.get_data_from_filename(filename)

            df = pd.read_csv(os.path.join(folders.predictions, filename))
            df['scale'] = scale
            df['exp_type'] = exp_type
            df['target_exp'] = target_exp
            df['itern'] = itern

            if '_avg' in filename:
                df.rename(columns={'avg': 'exp_value'}, inplace=True)
            elif '_rtx' in filename:
                df.rename(columns={'rtx': 'exp_value'}, inplace=True)
            else:
                assert False

            dfs.append(df)

        master_df = pd.concat(dfs, ignore_index=True)
        master_df.to_csv(os.path.join(
            folders.eqexp, 'predictions.csv'), index=False)
        return master_df

    def aggregate_model_results(data_df, scales=config.scales):
        ''' Calc and aggreagte models results over all simulated outcomes data '''

        assert isinstance(scales, list)
        assert isinstance(data_df, pd.core.frame.DataFrame)

        exp_types = config.exp_types
        target_exps = data_df.target_exp.unique()
        iterns = data_df.itern.unique()
        model_results = []

        for scale, exp_type, target_exp, itern in tqdmproduct(scales, exp_types, target_exps, iterns, desc='Calc models over simulations'):

            if itern>2:
                continue

            model = EqExpCore.get_model(
                data_df=data_df,
                scale=scale,
                exp_type=exp_type,
                target_exp=target_exp,
                itern=itern)

            p, ci_low, est, ci_high = EqExpCore.get_model_outputs(model=model)

            # Todo: slow code, create DF from dictionary to speed up
            model_results.append({
                'scale': scale,
                'exp_type': exp_type,
                'target_exp': target_exp,
                'itern': itern,
                'p': p,
                'ci_low': ci_low,
                'est': est,
                'ci_high': ci_high
            })

        stats_df = pd.DataFrame.from_dict(model_results)
        stats_df.to_csv(os.path.join(
            folders.eqexp, 'model_results.csv'), index=False)
        return stats_df

    def get_avg_across_iterns(df):
        ''' Calc avg p and effect estimate across iterations '''

        assert isinstance(df, pd.core.frame.DataFrame)

        scales = df.scale.unique()
        exp_types = df.exp_type.unique()
        target_exps = df.target_exp.unique()

        avg_df = copy.deepcopy(df)
        del avg_df['itern']
        del avg_df['p']
        del avg_df['est']
        del avg_df['ci_low']
        del avg_df['ci_high']
        avg_df['p'] = None
        avg_df['est'] = None
        avg_df['ci_low'] = None
        avg_df['ci_high'] = None
        avg_df.drop_duplicates(inplace=True)

        for scale, exp_type, target_exp in tqdmproduct(scales, exp_types, target_exps, desc='Calc avg model outcomes'):

            tmp = df.loc[
                (df.scale == scale) &
                (df.exp_type == exp_type) &
                (df.target_exp == target_exp)
            ]

            avg_df.loc[
                (df.scale == scale) &
                (df.exp_type == exp_type) &
                (df.target_exp == target_exp), 'p'] = round(tmp.p.mean(), 5)

            avg_df.loc[
                (df.scale == scale) &
                (df.exp_type == exp_type) &
                (df.target_exp == target_exp), 'est'] = round(tmp.est.mean(), 5)

            avg_df.loc[
                (df.scale == scale) &
                (df.exp_type == exp_type) &
                (df.target_exp == target_exp), 'ci_low'] = round(tmp.ci_low.mean(), 5)

            avg_df.loc[
                (df.scale == scale) &
                (df.exp_type == exp_type) &
                (df.target_exp == target_exp), 'ci_high'] = round(tmp.ci_high.mean(), 5)

        avg_df.to_csv(os.path.join(
            folders.eqexp, 'avg__model_results.csv'), index=False)

    def get_model(data_df, scale, exp_type, target_exp, itern):
        ''' Fits python model to simulated data '''

        assert isinstance(data_df, pd.core.frame.DataFrame)
        assert isinstance(scale, str)
        assert isinstance(exp_type, str)

        df = data_df.loc[
            (data_df.scale == scale) &
            (data_df.exp_type == exp_type) &
            (data_df.target_exp == target_exp) &
            (data_df.itern == itern) &
            (data_df.tp.isin(['wk0', 'wk6']))
        ]

        df = df.dropna()
        model = smf.mixedlm('score ~ trt*tp', data=df,
                            groups=df['patient_id']).fit()
        return model

    def get_model_outputs(model, key='trt[T.P]:tp[T.wk6]'):
        ''' Extract and returns specific model estimates '''

        p = round(model.pvalues[key], 6)
        ci_low = round(model.conf_int()[0][key], 3)
        ci_high = round(model.conf_int()[1][key], 3)
        est = round((model.conf_int()[0][key] + model.conf_int()[1][key])/2, 3)

        return p, ci_low, est, ci_high

    def get_data_from_filename(filename):
        ''' Extracts all info embedded into filenames '''

        assert isinstance(filename, str)

        # Get scale
        idx = filename.index('_scale')
        if filename[(idx+6):(idx+9)] == 'HAM':
            scale = 'HAMD'
        elif filename[(idx+6):(idx+9)] == 'BDI':
            scale = 'BDI'
        elif filename[(idx+6):(idx+9)] == 'MAD':
            scale = 'MADRS'
        elif filename[(idx+6):(idx+9)] == 'QID':
            scale = 'QIDS'
        elif filename[(idx+6):(idx+9)] == 'STA':
            scale = 'STAIT'
        elif filename[(idx+6):(idx+9)] == 'WEM':
            scale = 'WEMWBS'
        else:
            assert False

        # Get exp data
        if 'avg' in filename:
            assert 'rtx' not in filename
            exp_type = 'avg'
            idx = filename.index('avg')

        if 'rtx' in filename:
            assert 'avg' not in filename
            exp_type = 'rtx'
            idx = filename.index('rtx')

        target_exp = int(filename[(idx+3):(idx+6)])

        # Get iteration
        idx = filename.index('_iter')
        itern = int(filename[(idx+5):(idx+8)])
        return scale, exp_type, target_exp, itern


class EqExpFigures():

    def get_all_eqexp_figure(stats_df, avg_df, scales=config.scales, exp_types=config.exp_types, do_save=True):
        ''' Plots and saves all equal expectancy plots'''

        assert isinstance(stats_df, pd.core.frame.DataFrame)
        assert isinstance(avg_df, pd.core.frame.DataFrame)
        assert isinstance(scales, list)
        assert isinstance(exp_types, list)
        assert isinstance(do_save, bool)

        for scale, exp_type in tqdmproduct(scales, exp_types, desc='Plot results'):
            EqExpFigures.get_eqexp_figure(
                stats_df=stats_df,
                avg_df=avg_df,
                scale=scale,
                exp_type=exp_type,
                do_show=False,
                do_save=do_save)

    def get_eqexp_figure(stats_df, avg_df, scale='HAMD', exp_type='rtx', do_show=True, do_save=False):
        ''' Plots and saves a particular equal expectancy plots'''

        assert isinstance(stats_df, pd.core.frame.DataFrame)
        assert isinstance(avg_df, pd.core.frame.DataFrame)
        assert isinstance(exp_type, str)
        assert isinstance(do_show, bool)
        assert isinstance(do_save, bool)

        tmp = stats_df.loc[
            (stats_df.scale == scale) &
            (stats_df.exp_type == exp_type)
        ]

        tmp_avg = avg_df.loc[
            (avg_df.scale == scale) &
            (avg_df.exp_type == exp_type)
        ]

        fig = pyplt.figure(figsize=(6, 4))
        ax11 = fig.add_subplot(1, 1, 1)
        ax12 = ax11.twinx()
        ax12.grid(False)

        p_color = 'blue'
        est_color = 'green'

        ax11.set_ylabel('Treatment p-value', color=p_color,
                        fontdict=config.axislabel_fontdict)
        ax12.set_ylabel('Treatment estimate', color=est_color,
                        fontdict=config.axislabel_fontdict)
        ax11.set_title(scale, fontdict=config.title_fontdict)
        ax11.tick_params(axis='y', which='major', labelsize=14)
        ax12.tick_params(axis='y', which='major', labelsize=14)

        if exp_type == 'rtx':
            ax11.set_xlabel("Mean 'received treatment expectancy'",
                            fontdict=config.axislabel_fontdict)
        elif exp_type == 'avg':
            ax11.set_xlabel("Mean 'average trial expectancy'",
                            fontdict=config.axislabel_fontdict)
        else:
            assert False

        sns.lineplot(ax=ax11, y='p', x='target_exp', errorbar='ci',
                     n_boot=1024, data=tmp, color=p_color)
        sns.lineplot(ax=ax12, y='est', x='target_exp',
                     errorbar=None, data=tmp_avg, color=est_color)
        ax12.fill_between(
            x=tmp_avg.target_exp,
            y1=tmp_avg.ci_low,
            y2=tmp_avg.ci_high, color=est_color, alpha=0.35)

        if do_show:
            fig.show()
            input("Press any key to continue\n")

        if do_save:
            output_dir = os.path.abspath(os.path.join(folders.figures_export))
            output_fname = f'eqexp_{scale}_{exp_type}'
            dpi = 300

            fig.savefig(
                fname=os.path.join(output_dir, output_fname+'.png'),
                format='png',
                dpi=dpi,
            )

            fig.savefig(
                fname=os.path.join(output_dir, output_fname+'.svg'),
                format='svg',
                dpi=dpi,
            )
