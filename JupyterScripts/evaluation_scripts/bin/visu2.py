import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib as mpl
from evaluation_scripts.base import load_experiment, load_toml
from evaluation_scripts.visu import sort_legend_by
import seaborn as sns
import os
from itertools import product


def load_exp(exp):
    sns.set_style("darkgrid")
    _, _, _, pathes, Data, _ = load_experiment(
        exp, src="01b_Data")
    Data = Data[Data.exclude != True]
    color_dict = load_toml("phenotype_colors.toml")
    axs_dict = {"No treatment": (0, 0), "Mono B": (0, 1), "Mixing": (
        0, 2), "Mono A": (1, 0), "Cycling": (1, 1), "Combo": (1, 2)}
    return axs_dict, color_dict, pathes, Data


def load_model(name, pathes):
    model = pd.read_pickle(os.path.join(pathes["obj"], name)).reset_index().rename(
        columns={"r": "rep", "x_hat_sim": "phenotype"})
    return model


def plot_exp_data(ax, df, color_dict, scale=2):
    r = df.rep.unique()
    p = df.phenotype.unique()
    s = df.strategy.unique()[0]

    L = []
    for t in range(0, max(df.transfer_n+1)):
        if t not in df.transfer_n.unique():
            for phenotype, rep in product(p, r):
                L.append({"rep": rep, "phenotype": phenotype,
                         "f": None, "strategy": s, "transfer_n": t})
    fill = pd.DataFrame(L)
    df = pd.concat([fill, df], sort=True)
    sns.pointplot(data=df,
                  ax=ax,
                  x="transfer_n",
                  y="f",
                  scale=scale,
                  hue="phenotype",
                  errorbar=("pi", 100), join=False, dodge=0.2,
                  palette=color_dict,
                  )


def make_legend(rcparams, ax, fig):
    order = ["U", "S", "A_r", "B_r", "A&B", "AB_r", "Other",
             'exp data (95% pi)', 'validation sim', 'variation sim (95% pi)']
    ax.plot([], [], color=rcparams['lines.color'],
            linestyle="dashdot", label='validation sim', linewidth=3)
    handles, labels = ax.get_legend_handles_labels()

    h, l = make_exp_data_entry(
        'exp data (95% pi)', ax, rcparams['lines.color'])
    handles.append(h)
    labels.append(l)

    h, l = make_error_legend_entry(
        'variation sim (95% pi)', ax, rcparams['lines.color'])
    handles.append(h)
    labels.append(l)
    handles, labels = sort_legend_by(order, handles, labels)
    fig.legend(handles, labels, markerscale=6, bbox_to_anchor=(0.5, -0.1), loc='lower center', ncol=len(
        labels), prop={'size': 20}, frameon=False)


def make_legend_exp(rcparams, ax, fig, s):
    order = ["UI", "S", "A_r", "B_r", "A&B", "AB_r", "Fishy",
             'exp data (mean, min, max)']

    handles, labels = ax.get_legend_handles_labels()

    h, l = make_exp_data_entry(
        'exp data (mean, min, max)', ax, rcparams['lines.color'], s=3*s)
    handles.append(h)
    labels.append(l)

    handles, labels = sort_legend_by(order, handles, labels)
    fig.legend(handles, labels, markerscale=s, bbox_to_anchor=(0.5, 0), loc='lower center', ncol=len(
        labels), prop={'size': 20}, frameon=False)


def plot_exp_strategy(exp, rcparams, strategy, s=2):
    axs_dict, color_dict, pathes, df = load_exp(exp)
    cols = ["transfer_n", "strategy", "rep", "phenotype", "n"]
    df["n"] = True
    summary = df[cols].groupby(cols[:-1]).count().unstack(fill_value=0).stack()
    cols = ["transfer_n", "strategy", "rep", "n"]
    total = df[cols].groupby(cols[:-1]).count().unstack(fill_value=0).stack()
    for i, j in summary.iterrows():
        summary.loc[i, "f"] = summary.loc[i, "n"]/total.loc[i[:-1], "n"]
    summary = summary.reset_index()

    fontsize_small = 0.7*rcparams["font.size"]
    fontsize_large = 1.5*rcparams["font.size"]

    mpl.rcParams.update(**rcparams)
    fig, ax = plt.subplots(1, 1, figsize=(
        25, 15))

    plot_exp_data(
        ax, summary[summary.strategy.isin([strategy])], color_dict, scale=s)

    ax.get_legend().remove()
    fig.suptitle('TS-Plasmid: '+exp, fontsize=fontsize_large)
    ax.tick_params(axis='both', which='major', labelsize=fontsize_small)
    make_legend_exp(rcparams, ax, fig, s=s)
    ax.set_xlabel('transfer', fontsize=fontsize_small)
    ax.set_ylabel(
        'fraction of patients with phenotype', fontsize=fontsize_small)

    plt.savefig(os.path.join(
        pathes["analysis"], exp+"_exp_data_plot_"+strategy+".pdf"), bbox_inches='tight', dpi=300)


def make_exp_data_entry(entry, ax, color, s=2):
    p1, = ax.plot([], [], c=color, marker="|", markersize=5*s,
                  linestyle='None')  # notice the comma!
    p2, = ax.plot([], [], c=color, marker="o", markersize=s,
                  linestyle='None')  # notice the comma!
    handle = ((p1, p2),)
    label = entry
    return handle, label


def lineplot(ax, df, color_dict, s=100, errorbar=None, line_style="-"):
    l0 = len(ax.lines)
    sns.lineplot(data=df, ax=ax,
                 x="transfer_n",
                 y="f",
                 hue="phenotype",
                 palette=color_dict,
                 errorbar=errorbar,
                 legend=None,
                 )
    for l in ax.lines[l0:]:
        l.set_linestyle(line_style)


def set_rcparams(style="white", fontsize=30, box=False):
    black = plt.style.library["dark_background"].copy()
    white = black.copy()
    for key in black:
        if black[key] == "black":
            white.update({key: "white"})
        elif black[key] == "white":
            white.update({key: "black"})
    if style == "white":
        ek = "black"
        rcparams = white
    else:
        rcparams = black
        ek = "white"

    rcparams.update({**rcparams, 'font.size': fontsize})
    if not box:
        rcparams['axes.edgecolor'] = rcparams['axes.facecolor']

    return removekey(rcparams, "grid.color")


def removekey(d, key):
    r = dict(d)
    del r[key]
    return r


def sort_legend_by(order, handles, labels):
    rank = list(range(len(order)))
    order_dict = dict(zip(order, rank))
    tmp = pd.DataFrame()
    tmp["handles"] = handles
    tmp["labels"] = labels
    tmp["rank"] = tmp.apply(lambda x: order_dict[x["labels"]], axis=1)
    tmp.sort_values("rank", inplace=True)
    return list(tmp.handles), list(tmp.labels)


def make_error_legend_entry(entry, ax, color):
    if color == "black":
        alpha = 0.2
    else:
        alpha = 0.4
    p1, = ax.plot([], [], c=color, linewidth=3,
                  linestyle="solid")  # notice the comma!
    p2, = ax.plot([], [], c=color, linestyle="-", linewidth=20,
                  alpha=alpha)  # notice the comma!
    handle = ((p1, p2),)
    label = entry
    return handle, label


def plot_exp(exp, rcparams, n=1, name_add="", filter="_unfiltered", grid_alpha=0.2, x_tick_space=2):
    fontsize_small = 0.7*rcparams["font.size"]
    fontsize_large = 1.5*rcparams["font.size"]
    args = get_plot_data(exp, filter, x_tick_space)
    axs_dict = args[0]

    mpl.rcParams.update(**rcparams)
    fig, axs = plt.subplots(2, 3, figsize=(
        25, 15), sharex=True, sharey=True, layout="constrained")
    fig.suptitle('TS-Plasmid: '+exp, fontsize=rcparams["font.size"])

    for strategy in axs_dict.keys():

        j, k = axs_dict[strategy]
        ax = axs[axs_dict[strategy]]

        plot_strategy(ax, args, rcparams, grid_alpha, n, strategy)

        if j + k == 0:
            make_legend(rcparams, ax, fig)

        if j == 1:
            ax.set_xlabel('transfer', fontsize=fontsize_small)

        if k == 0:
            ax.set_ylabel(
                'fraction of patients with phenotype [%]', fontsize=fontsize_small)

        if k != 0:
            ax.set_ylabel("")
        else:
            ax.set_ylabel("phenotype frequency")

        if j == 0:
            ax.set_xlabel("")
        else:
            ax.set_xlabel("transfer")

    fig.get_layout_engine().set(hspace=0.1)
    pathes = args[2]
    plt.savefig(os.path.join(
        pathes["analysis"], exp+name_add+"_plot.pdf"), bbox_inches='tight', dpi=300)


def plot_strategy(ax, args, rcparams, grid_alpha, n, strategy):
    axs_dict, color_dict, pathes, summary, model_var, model_val, x_ticks = args
    plot_exp_data(
        ax, summary[summary.strategy.isin([strategy])], color_dict)

    lineplot(ax, model_val[model_val.strategy.isin(
        [strategy])], color_dict, line_style="-.")

    for _ in range(n):
        lineplot(ax, model_var[model_var.strategy.isin(
            [strategy])], color_dict, errorbar="pi")

    ax.get_legend().remove()
    ax.grid(True, color=rcparams['xtick.color'], alpha=grid_alpha)
    ax.set_xticks(x_ticks)
    ax.set_title(strategy, fontsize=rcparams["font.size"])

    for i, tick in enumerate(ax.xaxis.get_major_ticks()):
        if i % 2 == 0:
            tick.label1.set_fontsize(18)
        else:
            tick.label1.set_visible(False)

    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(18)


def get_plot_data(exp, filter, x_tick_space):
    axs_dict, color_dict, pathes, df = load_exp(exp)
    cols = ["transfer_n", "strategy", "rep", "phenotype", "n"]
    df["n"] = True
    summary = df[cols].groupby(cols[:-1]).count().unstack(fill_value=0).stack()
    cols = ["transfer_n", "strategy", "rep", "n"]
    total = df[cols].groupby(cols[:-1]).count().unstack(fill_value=0).stack()
    for i, j in summary.iterrows():
        summary.loc[i, "f"] = summary.loc[i, "n"]/total.loc[i[:-1], "n"]
    summary = summary.reset_index()
    model_val = load_model("05c_bbb_sim_val_av"+filter+".pkl", pathes)
    model_var = load_model("05b_bbb_sim_altern_av"+filter+".pkl", pathes)
    T = df.transfer_n.unique()
    x_ticks = T[::x_tick_space]
    return axs_dict, color_dict, pathes, summary, model_var, model_val, x_ticks


def plot_single(exp, rcparams, strategy, n=1, name_add="", filter="_unfiltered", grid_alpha=0.2, x_tick_space=2):
    fontsize_small = 0.7*rcparams["font.size"]
    fontsize_large = 1.5*rcparams["font.size"]
    args = get_plot_data(exp, filter, x_tick_space)
    pathes = args[2]
    mpl.rcParams.update(**rcparams)
    fig, ax = plt.subplots(1, 1, figsize=(
        25, 15))
    fig.suptitle('TS-Plasmid: '+exp, fontsize=fontsize_large)
    plot_strategy(ax, args, rcparams, grid_alpha, n, strategy)
    ax.tick_params(axis='both', which='major', labelsize=fontsize_small)
    make_legend(rcparams, ax, fig)
    ax.set_xlabel('transfer', fontsize=fontsize_small)
    ax.set_ylabel(
        'fraction of patients with phenotype', fontsize=fontsize_small)

    plt.savefig(os.path.join(
        pathes["analysis"], exp+name_add+"_plot_"+strategy+".pdf"), bbox_inches='tight', dpi=300)
