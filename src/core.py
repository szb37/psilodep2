"""
:Author: Balazs Szigeti {szb37 AT pm DOT me}
:Copyright: 2022, DrugNerdsLab
:License: MIT
"""

import matplotlib.pyplot as pyplt
import src.folders as folders
import src.config as fig_config
import seaborn as sns
import pandas as pd
import os


class Plots():

    def boxplots_per_treatment(df, variable, y_label, show=False, save=True):
        ''' Boxplots of baseline between-arms differences'''

        assert isinstance(df, pd.core.frame.DataFrame)
        assert isinstance(variable, str)
        assert variable in df.columns
        assert isinstance(y_label, str)
        assert isinstance(show, bool)
        assert isinstance(save, bool)

        # Format data
        df = df[['patient_id', 'trt', variable]]
        df.drop_duplicates(inplace=True)
        df = df[~df.isnull().any(axis=1)]

        # Setup figure
        fig = pyplt.figure(figsize=(4, 4), dpi=200)
        ax = fig.add_subplot(1, 1, 1)

        sns.boxplot(ax=ax, data=df, x='trt', y=variable, order=['E', 'P'])

        ax.set_xlabel(None)
        ax.set_xticklabels(['Escitalopram arm', 'Psilocybin arm'],
                           fontdict=fig_config.ticklabel_fontdict)
        ax.set_ylabel(y_label, fontdict=fig_config.axislabel_fontdict)
        ax.tick_params(axis='both', which='major',
                       labelsize=fig_config.tick_labelsize)

        Helpers.save_figure(
            fig=fig,
            show=show,
            save=save,
            out_folder=os.path.join(folders.figures_export),
            fname=f'boxplot_{variable}',
        )

    def regression(df, variable, x_label, xlim=None, show=False, save=True):
        ''' Regression line of baseline variable vs outcome'''

        assert isinstance(df, pd.core.frame.DataFrame)
        assert isinstance(variable, str)
        assert variable in df.columns
        assert isinstance(x_label, str)
        assert isinstance(show, bool)
        assert isinstance(save, bool)
        assert (xlim is None or isinstance(xlim, list))

        # Format data
        df = df.loc[(df.tp == 'wk6')]
        df = df[['patient_id', 'trt', 'scale', 'score_chg', variable]]
        df.drop_duplicates(inplace=True)
        df = df[~df.isnull().any(axis=1)]

        scatter_kws = {'alpha': 0.5, 's': 7, 'edgecolors': None}
        robust = True
        scatter = True
        truncate = False

        for scale in ['hamd', 'bdi', 'madrs', 'qids', 'stait', 'wemwbs']:

            fig = pyplt.figure(figsize=(5, 4), dpi=200)
            ax = fig.add_subplot(1, 1, 1)
            if xlim is not None:
                ax.set_xlim(xlim)
            ax.set_title(scale.upper(), fontdict=fig_config.title_fontdict)
            ax.tick_params(axis='both', which='major',
                           labelsize=fig_config.tick_labelsize)

            df_esc = df.loc[(df.trt == 'E') & (df.scale == scale)]
            df_psi = df.loc[(df.trt == 'P') & (df.scale == scale)]

            sns.regplot(
                ax=ax, x_ci='ci', data=df_esc, x=variable, y='score_chg', color='blue',
                scatter=scatter, scatter_kws=scatter_kws, robust=robust, truncate=truncate,
                line_kws={'label': 'Escitalopram arm'}
            )
            sns.regplot(
                ax=ax, x_ci='ci', data=df_psi, x=variable, y='score_chg', color='red',
                scatter=scatter, scatter_kws=scatter_kws, robust=robust, truncate=truncate,
                line_kws={'label': 'Psilocybin arm'}
            )

            ax.set_xlabel(x_label, fontdict=fig_config.axislabel_fontdict)
            ax.set_ylabel('Δ score (wk6-wk0)',
                          fontdict=fig_config.axislabel_fontdict)

            ax.legend(fontsize=fig_config.legend_fontsize)

            Helpers.save_figure(
                fig=fig,
                show=show,
                save=save,
                out_folder=os.path.join(folders.figures_export),
                fname=f'regression_{scale}_{variable}',
            )


class Helpers():

    def save_figure(fig, show, save, out_folder, fname):
        ''' Save figure '''

        assert isinstance(show, bool)
        assert isinstance(save, bool)
        assert isinstance(out_folder, str)
        assert isinstance(fname, str)

        if show:
            fig.show()
            input("Press any key to continue\n")

        if save:
            dpi = 300

            fig.savefig(
                fname=os.path.join(out_folder, fname+'.png'),
                format='png',
                dpi=dpi,
            )

            fig.savefig(
                fname=os.path.join(out_folder, fname+'.svg'),
                format='svg',
                dpi=dpi,
            )
