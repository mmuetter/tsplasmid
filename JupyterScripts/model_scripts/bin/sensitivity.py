from evaluation_scripts.base import get_base_pathes, load_toml
import os
import pandas as pd
from matplotlib import pyplot as plt
from itertools import product
from matplotlib.ticker import MaxNLocator


def plot_results(summary, strategies, colors, handles, style="white", lw=5,  g=None, exp=False, case="winner", x="R0", y="single_res", x_label="$R_0$", y_label="$c_A + c_B$", marker="o"):
    black = plt.style.library["dark_background"].copy()
    white = black.copy()
    for key in black:
        if black[key] == "black":
            white.update({key: "white"})
        elif black[key] == "white":
            white.update({key: "black"})
    if style == "white":
        ek = "black"
        plt.rcParams.update({**white, 'font.size': 50})
    else:
        plt.rcParams.update({**black, 'font.size': 50})
        ek = "white"

    if g:
        fig, ax = g
    else:
        fig, ax = plt.subplots(figsize=(30, 20))

    if exp == False:
        s = "s"
        alphas = summary.alpha.unique()

        for strategy, alpha in product(strategies, alphas):
            df = summary[(summary[case] == strategy)
                         & (summary.alpha == alpha)]

            g = ax.scatter(x=x,
                           y=y,
                           s=s,
                           color=colors[strategy],
                           marker=marker,
                           alpha=alpha,
                           data=df)
        ax.legend(labels=strategies, loc=(1.04, 0))
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.legend(handles=handles,
                  markerscale=2,
                  bbox_to_anchor=(0.5, -0.15),
                  loc='lower center',
                  ncol=6,
                  prop={'size': 20},
                  frameon=False)
    else:
        if case == "winner":
            f = 10
        else:
            f = 1
        for _, df in summary.iterrows():
            g = ax.scatter(x=x,
                           y=y,
                           s=df["s"],
                           color=colors[df[case]],
                           marker="d",
                           alpha=df["alpha"],
                           data=df,
                           edgecolors=ek)

            X = df[x]
            Y = df[y]
            ax.annotate(df["name"], (X+f*0.02, Y), color=ek)

    return fig, ax


def plot_results_2(summary, strategies, colors, handles, style="white", lw=5,  g=None, exp=False, case="winner", x="R0", y="single_res", x_label="$R_0$", y_label="$c_A + c_B$", marker="o", s=1000, x_lim=(0, 10), y_lim=(0, 1)):
    black = plt.style.library["dark_background"].copy()
    white = black.copy()
    for key in black:
        if black[key] == "black":
            white.update({key: "white"})
        elif black[key] == "white":
            white.update({key: "black"})
    if style == "white":
        ek = "black"
        plt.rcParams.update({**white, 'font.size': 50})
    else:
        plt.rcParams.update({**black, 'font.size': 50})
        ek = "white"

    if g:
        fig, ax = g
    else:
        fig, ax = plt.subplots(figsize=(30, 20))

    if exp == False:
        alphas = summary.alpha.unique()

        for strategy, alpha in product(strategies, alphas):
            df = summary[(summary[case] == strategy)
                         & (summary.alpha == alpha)]

            g = ax.scatter(x=x,
                           y=y,
                           s=s,
                           color=colors[strategy],
                           marker=marker,
                           alpha=alpha,
                           data=df)
        ax.legend(labels=strategies, loc=(1.04, 0))
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.legend(handles=handles,
                  markerscale=2,
                  bbox_to_anchor=(0.5, -0.2),
                  loc='lower center',
                  ncol=6,
                  prop={'size': 40},
                  frameon=False)
    else:
        if case == "winner":
            f = 10
        else:
            f = 1
        for _, df in summary.iterrows():
            g = ax.scatter(x=x,
                           y=y,
                           s=s,
                           color=colors[df[case]],
                           marker=marker,
                           alpha=1,
                           data=df,
                           edgecolors=ek)

            X = df[x]
            Y = df[y]
            ax.annotate(df["name"], (X+f*0.02, Y-f*0.005), color=ek)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=6, prune='lower'))

    ax.set_box_aspect(1)
    return fig, ax


def get_summary(file_name, case, crit):
    _, base = get_base_pathes()
    path = os.path.join(base, "sensitivity_sim")
    df = pd.read_pickle(os.path.join(path, file_name))
    colors = load_toml("strategy_colors.toml")
    summary = df[df[case].isnull() == False].copy()
    summary = summary[(summary.margin.isnull() == False) &
                      (summary.mean_margin.isnull() == False)]
    summary["single_res"] = summary.c_A + summary.c_B
    summary["R0"] = summary.infection/summary.turnover
    summary["alpha"] = round(summary[crit] / max(summary[crit]), 2)

    strategies = []
    handles = []
    for strategy in summary[case].unique():
        strategies.append(strategy)
        handles.append(plt.Line2D([], [], color=colors[strategy],
                                  marker="o", markersize=25, linewidth=0, label=strategy))
    return summary, handles, strategies, colors, path


def make_triangle(case):
    A = [-1, 1]
    B = [1, 1]
    C = [-1, -1]
    D = [1, -1]
    if case == "upper":
        marker = [A, B, C]
    elif case == "lower":
        marker = [C, B, D]
    return marker


def prep_data(file_name, s0=5000):
    _, base = get_base_pathes()
    path = os.path.join(base, "sensitivity_sim")
    summary = pd.read_pickle(os.path.join(path, file_name))
    colors = load_toml("strategy_colors.toml")
    summary_l, handles_l, strategies_l = prep_frames(
        summary, colors, s0, "loser")
    summary, handles, strategies = prep_frames(summary, colors, s0)
    return strategies, strategies_l, handles, handles_l, summary, summary_l, colors, path


def prep_frames(df, colors, s0, case="winner"):
    if case == "winner":
        add = ""
    else:
        add = "_l"
    summary = df[df[case].isnull() == False].copy()
    summary = summary[(summary.margin.isnull() == False) &
                      (summary.mean_margin.isnull() == False)]
    summary["single_res"] = summary.A_r + summary.B_r
    summary["R0"] = summary.infection/summary.turnover
    summary["s"] = s0 * summary["mean_margin"+add] / \
        max(summary["mean_margin"+add])
    summary["alpha"] = round((summary["margin"+add]) /
                             max(summary["margin"+add]), 2)
    strategies = []
    handles = []
    for strategy in summary["winner"].unique():
        strategies.append(strategy)
        handles.append(plt.Line2D([], [], color=colors[strategy],
                                  marker="o", markersize=25, linewidth=0, label=strategy))
    return summary, handles, strategies
