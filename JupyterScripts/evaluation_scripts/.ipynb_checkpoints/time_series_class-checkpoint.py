from evaluation_scripts.legend_class import Legend
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os


class TimeSeries:
    def __init__(self, experiment):
        self.experiment = experiment
        self.order = self.experiment.strategies
        self.pathes = experiment.pathes
        self.summarize_data()
        self.model_data = []
        self.model_names = []
        self.model_dict = {}
        self.T = experiment.T

    def summarize_data(self):
        df = self.experiment.data.copy()
        cols = ["transfer_n", "strategy", "rep", "phenotype", "n"]
        df["n"] = True
        summary = df[cols].groupby(
            cols[:-1]).count().unstack(fill_value=0).stack()
        cols = ["transfer_n", "strategy", "rep", "n"]
        total = df[cols].groupby(
            cols[:-1]).count().unstack(fill_value=0).stack()
        for i, j in summary.iterrows():
            summary.loc[i, "f"] = summary.loc[i, "n"]/total.loc[i[:-1], "n"]
        self.summary = summary.reset_index()

    def plot_strategy(self, strategy, aes=None, ax=None, sans_poisson=True):
        if sans_poisson:
            df = self.summary[self.summary.phenotype != "Other"].copy()
        else:
            df = self.summary.copy()
        self.strat_data = strat_data = df[df.strategy.isin([strategy])]
        
        if ax is None:
            _, ax = plt.subplots(1)
    
        # Experimental Data
        self.strategy = strategy
        sns.pointplot(
            ax=ax,
            data=strat_data,
            hue="phenotype",
            y="f",
            x="transfer_n",
            order=self.T,
            errorbar=("pi", 100),
            join=False,
            dodge=0.2,
            palette=aes.phenotype_colors if aes else None  # Check for aes
        )
    
        # Models
        for model, model_case in zip(self.model_data, self.model_dict.values()):
            plot_model(model, strategy, model_case, ax, aes)
    
        # Aesthetics
        if aes:
            mpl.rcParams.update(**aes.rcparams)
            ax.set_title(strategy, fontsize=aes.fontsize_large)
            ax.get_legend().remove()
            x_ticks = self.experiment.T[::aes.x_tick_space]
            ax.grid(
                aes.grid, color=aes.rcparams['xtick.color'], alpha=aes.grid_alpha)
            ax.set_xticks(x_ticks)
            
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(aes.fontsize_small)
        else:
            ax.set_title(strategy)
            ax.grid(True)


    def load_models(self, model_dict, model_new=True):
        self.model_dict = model_dict
        self.model_path = self.pathes["simulations"] if model_new else self.pathes["obj"]
        self.model_new = model_new
        self.model_names = self.model_dict.keys()
        for model_name in self.model_names:
            model = self.load_model(model_name).rename(
                columns={"x_hat": "phenotype"})
            self.model_data.append(model)

    def plot_all(self, aes, sans_poisson=True, panel_style = False):
        experiment = self.experiment

        fig, axs = plt.subplots(2, 3, figsize=(
            25, 15), sharex=True, sharey=True, layout="constrained")

        if aes.title:
            fig.suptitle('TS-Plasmid: '+experiment.exp,
                         fontsize=aes.fontsize_large)
    
        subplot_labels = ['(A)', '(B)', '(C)', '(D)', '(E)', '(F)']  # Labels for the subplots
        label_index = 0

        for strategy in experiment.axs_dict.keys():

            j, k = experiment.axs_dict[strategy]
            ax = axs[experiment.axs_dict[strategy]]

            self.plot_strategy(strategy, aes, ax=ax, sans_poisson=sans_poisson)
            if panel_style:
                ax.text(.1, 1.1, subplot_labels[label_index], transform=ax.transAxes, 
                        fontsize=aes.fontsize_large, va='top', ha='right')
                label_index += 1

            if j + k == 0:
                legend = Legend(fig, ax, aes)
                if "exp" in aes.legend.keys():
                    legend.make_exp_data_entry()
                if "var" in aes.legend.keys():
                    legend.make_var_legend()
                if "val" in aes.legend.keys():
                    legend.make_val_legend()
                legend.show()

            if j == 1:
                ax.set_xlabel(aes.xlabel, fontsize=aes.fontsize_small)

            if k == 0:
                ax.set_ylabel(
                    aes.ylabel, fontsize=aes.fontsize_small)

            if k != 0:
                ax.set_ylabel("")
            else:
                ax.set_ylabel(aes.ylabel)

            if j == 0:
                ax.set_xlabel("")
            else:
                ax.set_xlabel(aes.xlabel)

        fig.get_layout_engine().set(hspace=0.1)

    def load_model(self, name):
        model = pd.read_pickle(os.path.join(
            self.model_path, name)).reset_index()
        if self.model_new:
            model.rename(
                columns={"r": "rep", "x_hat": "phenotype"})
        else:
            model.rename(
                columns={"r": "rep", "x_hat_sim": "phenotype"})
        return model

    def save_figure(self, name, path="analysis"):
        fig_path = os.path.join(self.pathes[path], name)
        plt.savefig(fig_path, bbox_inches='tight')


def plot_model(model, strategy, case, ax, aes):
    if case == "var":
        for _ in range(aes.plot_n_times):
            lineplot(ax, model[model.strategy.isin([strategy])],
                     aes.phenotype_colors, errorbar="pi")
    elif case == "val":
        lineplot(ax, model[model.strategy.isin([strategy])],
                 aes.phenotype_colors, line_style="-.")
    else:
        print("case has to be var or val")


def lineplot(ax, df, color_dict, s=100, errorbar=None, line_style="-"):
    l0 = len(ax.lines)
    sns.lineplot(data=df, ax=ax,
                 x="transfer_n",
                 y="f",
                 hue="phenotype",
                 palette=color_dict,
                 errorbar=errorbar,
                 legend=None)
    for l in ax.lines[l0:]:
        l.set_linestyle(line_style)
