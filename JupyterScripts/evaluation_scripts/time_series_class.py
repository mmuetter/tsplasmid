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
        summary = df[cols].groupby(cols[:-1]).count().unstack(fill_value=0).stack()
        cols = ["transfer_n", "strategy", "rep", "n"]
        total = df[cols].groupby(cols[:-1]).count().unstack(fill_value=0).stack()
        for i, j in summary.iterrows():
            summary.loc[i, "f"] = summary.loc[i, "n"] / total.loc[i[:-1], "n"]
        self.summary = summary.reset_index()

    def strategy_plotter(self, strat_data, ax=None, palette=None):
        if not ax:
            _, ax = plt.subplots()

        sns.pointplot(
            ax=ax,
            data=strat_data,
            hue="phenotype",
            y="f",
            x="transfer_n",
            order=self.T,
            errorbar=("pi", 100),
            linestyle="none",  # Replace join=False with linestyle='none'
            dodge=0.2,
            palette=palette,
        )

        if ax.legend_:
            ax.legend_.remove()

    def plot_strategy(self, strategy, aes=None, ax=None, sans_poisson=True):
        if sans_poisson:
            df = self.summary[self.summary.phenotype != "Other"].copy()
        else:
            df = self.summary.copy()
        strat_data = df[df.strategy.isin([strategy])]

        if ax is None:
            _, ax = plt.subplots(1)

        # Experimental Data
        self.strategy_plotter(
            strat_data,
            ax,
            palette=(
                aes.phenotype_colors
                if aes and hasattr(aes, "phenotype_colors")
                else None
            ),
        )

        # Models
        self.plot_models(strategy, ax, aes.phenotype_colors, aes=aes)

        # Aesthetics
        if aes and hasattr(aes, "rcparams"):
            mpl.rcParams.update(**aes.rcparams)
        ax.set_title(
            strategy,
            fontsize=(
                aes.fontsize_large if aes and hasattr(aes, "fontsize_large") else 12
            ),
        )

        if aes and hasattr(aes, "x_tick_space"):
            x_ticks = self.experiment.T[:: aes.x_tick_space]
            ax.set_xticks(x_ticks)
        if aes and hasattr(aes, "grid"):
            ax.grid(aes.grid, color=aes.rcparams["xtick.color"], alpha=aes.grid_alpha)
        else:
            ax.grid(True)

        x_min, x_max = strat_data["transfer_n"].min(), strat_data["transfer_n"].max()
        y_min, y_max = strat_data["f"].min(), strat_data["f"].max()

        # Adjust limits with a small margin
        ax.set_xlim(x_min - 0.025 * (x_max - x_min), x_max + 0.025 * (x_max - x_min))
        ax.set_ylim(y_min - 0.05 * (y_max - y_min), y_max + 0.05 * (y_max - y_min))

        if aes and hasattr(aes, "fontsize"):
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(aes.fontsize)
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(aes.fontsize)

    def plot_models(self, strategy, ax, palette, aes=None):
        for model, model_case in zip(self.model_data, self.model_dict.values()):
            plot_model(model, strategy, model_case, ax, palette, aes)

    def load_models(self, model_dict, model_new=True):
        self.model_dict = model_dict
        self.model_path = (
            self.pathes["simulations"] if model_new else self.pathes["obj"]
        )
        self.model_new = model_new
        self.model_names = self.model_dict.keys()
        for model_name in self.model_names:
            model = self.load_model(model_name).rename(columns={"x_hat": "phenotype"})
            self.model_data.append(model)

    def plot_all(
        self,
        aes=None,
        sans_poisson=True,
        panel_style=False,
        layout=(2, 3),
        figsize=(25, 15),
        x_text=0.1,
        y_text=1.08,
    ):
        experiment = self.experiment

        fig, axs = plt.subplots(
            *layout, figsize=figsize, sharex=True, sharey=True, layout="constrained"
        )

        if aes and hasattr(aes, "title"):
            fig.suptitle(
                "TS-Plasmid: " + experiment.exp,
                fontsize=aes.fontsize_large if hasattr(aes, "fontsize_large") else 16,
            )

        subplot_labels = ["A)", "B)", "C)", "D)", "E)", "F)"]
        label_index = 0

        for strategy in experiment.axs_dict.keys():
            j, k = experiment.axs_dict[strategy]
            ax = axs[experiment.axs_dict[strategy]]

            self.plot_strategy(strategy, aes, ax=ax, sans_poisson=sans_poisson)

            if panel_style:
                ax.text(
                    x_text,
                    y_text,
                    subplot_labels[label_index],
                    transform=ax.transAxes,
                    fontsize=(
                        aes.fontsize_large
                        if aes and hasattr(aes, "fontsize_large")
                        else 12
                    ),
                    va="top",
                    ha="right",
                    weight="bold",
                )
                label_index += 1

            # First legend (primary entries like phenotype)
            if j + k == 0 and aes:
                primary_legend = Legend(fig, ax, aes)
                primary_legend.show()  # Show the first legend

            # Second legend (optional entries like exp, val, var)
            if j + k == 0 and aes:
                secondary_legend = Legend(fig, ax, aes, start_empty=True)
                # Optionally create entries for exp, val, var and call them second_row
                if hasattr(aes, "legend"):
                    if "exp" in aes.legend.keys():
                        secondary_legend.make_exp_data_entry()
                    if "var" in aes.legend.keys():
                        secondary_legend.make_var_legend()
                    if "val" in aes.legend.keys():
                        secondary_legend.make_val_legend()
                secondary_legend.show(i=1)

            ax.set_xlabel(
                aes.xlabel if aes and hasattr(aes, "xlabel") else "X-axis",
                fontsize=(aes.fontsize if aes and hasattr(aes, "fontsize") else 10),
            )
            ax.set_ylabel(
                aes.ylabel if aes and hasattr(aes, "ylabel") else "Y-axis",
                fontsize=(aes.fontsize if aes and hasattr(aes, "fontsize") else 10),
            )

        fig.get_layout_engine().set(hspace=0.1)
        return fig, axs

    def load_model(self, name):
        model = pd.read_pickle(os.path.join(self.model_path, name)).reset_index()
        if self.model_new:
            model.rename(columns={"r": "rep", "x_hat": "phenotype"})
        else:
            model.rename(columns={"r": "rep", "x_hat_sim": "phenotype"})
        return model

    def save_figure(self, name, path="analysis"):
        fig_path = os.path.join(self.pathes[path], name)
        plt.savefig(fig_path, bbox_inches="tight", dpi=300)


def plot_model(model, strategy, case, ax, palette, aes=None):
    if case == "var":
        for _ in range(aes.plot_n_times if aes and hasattr(aes, "plot_n_times") else 1):
            lineplot(
                ax,
                model[model.strategy.isin([strategy])],
                color_dict=palette,
                errorbar="pi",
            )
    elif case == "val":
        lineplot(
            ax,
            model[model.strategy.isin([strategy])],
            color_dict=palette,
            line_style="-.",
        )
    else:
        print("case has to be var or val")


def lineplot(ax, df, color_dict=None, s=100, errorbar=None, line_style="-"):
    l0 = len(ax.lines)
    sns.lineplot(
        data=df,
        ax=ax,
        x="transfer_n",
        y="f",
        hue="phenotype",
        palette=color_dict,
        errorbar=errorbar,
        legend=None,
    )
    for l in ax.lines[l0:]:
        l.set_linestyle(line_style)
