import seaborn as sns
import matplotlib.pyplot as plt
from evaluation_scripts.experiment_class import Experiment


class Heatmap:
    def __init__(self,
                 figsize=(20, 12),
                 cmap=sns.diverging_palette(366, 250,  as_cmap=True),
                 xticks_step=3,
                 yticks_step=3,
                 xlabel="$R_0$",
                 ylabel="$c_A + c_B$",
                 cbar_label="advantage",
                 fontsize_large=24,
                 fontsize=18):
        self.figsize = figsize
        self.cmap = cmap
        self.xticks_step = xticks_step
        self.yticks_step = yticks_step
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.cbar_label = cbar_label
        self.fontsize_large = fontsize_large
        self.fontsize = fontsize

    def plot_heatmap(self, df, dates, labels, show_colorbar=True):
        heatmap = df.pivot(index="res", columns="R0", values="score_max")
        heatmap = heatmap.sort_index(ascending=False)
        plt.figure(figsize=self.figsize)
        self.g = g = sns.heatmap(heatmap, cmap=self.cmap, vmin=-0.5, vmax=0.5,
                                 square=True, cbar=show_colorbar, cbar_kws={'label': self.cbar_label})

        # Set axis label font size and labels
        plt.xlabel(self.xlabel, fontsize=self.fontsize_large)
        plt.ylabel(self.ylabel, fontsize=self.fontsize_large)

        plt.xticks(fontsize=self.fontsize)
        plt.yticks(fontsize=self.fontsize)

        # Show every xticks_step tick label for the x-axis, starting from the first one
        for index, label in enumerate(g.get_xticklabels()):
            if index % self.xticks_step != 0:
                label.set_visible(False)

        # Show every yticks_step tick label for the y-axis, starting from the first one
        for index, label in enumerate(g.get_yticklabels()):
            if index % self.yticks_step != 0:
                label.set_visible(False)

        # Increase the size of the colorbar ticks and labels
        if show_colorbar:
            cbar = g.collections[0].colorbar
            cbar.ax.tick_params(labelsize=14)
            cbar.ax.set_ylabel(self.cbar_label, fontsize=16)

        for date, label in zip(dates, labels):
            exp = get_exp_stats(date)
            self.plt_exp(exp, label)

    def plt_exp(self, exp_dict, name, average=False):
        R0 = exp_dict["R0"]
        res = exp_dict["res"]
        y_ticks = self.g.axes.yaxis.get_ticklabels()
        y = get_tick_position(y_ticks, "y", res)

        x_ticks = self.g.axes.xaxis.get_ticklabels()
        x = get_tick_position(x_ticks, "x", R0)

        if average:
            color = self.cmap(exp_dict["score_av"] + 0.5)
        else:
            color = self.cmap(exp_dict["score"] + 0.5)
        plt.plot(x, y, "D", markersize=16,
                 markeredgecolor="black", color=color)

        plt.text(x+0.5, y+0.5, name, fontsize=self.fontsize)


def get_exp_stats(date):
    exp = Experiment(date)
    rates = exp.exp_pars["rates"][0]
    R0 = rates["infection"] / rates["turnover"]
    res = exp.community_distribution["A_r"] + exp.community_distribution["B_r"]
    T = max(exp.T)
    exp.make_time_series()
    data = exp.time_series.summary
    data = data[data.transfer_n.isin(
        [T-3, T-2, T-1, T]) & data.phenotype.isin(["U"])]
    cols = ["strategy", "f"]
    means = data[cols].groupby(cols[:-1]).mean()
    comparison = exp.strategies
    comparison.remove("Combo")
    best_comp = max(means.loc[comparison, "f"])
    score = means.loc["Combo", "f"] - best_comp
    score_av = means.loc["Combo", "f"] - means.loc[comparison, "f"].mean()
    return {"R0": R0,
            "res": res,
            "score": score,
            "score_av": score_av
            }


def interpolate(X, Y, x):
    i1 = X.index(min(X))
    i2 = X.index(max(X))
    x1 = X[i1]
    x2 = X[i2]
    y1 = Y[i1]
    y2 = Y[i2]
    m = (y2-y1)/(x2 - x1)
    y = y1 + m*x
    if y < 0:
        y = 0.4
    return y


def get_tick_position(ticks, axes, value):
    if axes == "y":
        entry = 1
    elif axes == "x":
        entry = 0

    labels = [float(ticks[0].get_text()), float(ticks[-1].get_text())]
    positions = [ticks[0].get_position()[entry],
                 ticks[-1].get_position()[entry]]
    return interpolate(labels, positions, value)
