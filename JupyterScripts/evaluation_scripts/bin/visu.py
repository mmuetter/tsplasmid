import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import random


def get_min_max_for_phenotypes(df, num_of_patients, phenotypes):
    f_min = {}
    f_max = {}
    f_mean = {}
    T = df.transfer_n.unique()
    T.sort()
    for phenotype in phenotypes:
        min_v = []
        max_v = []
        mean_v = []

        for t in T:
            # for r in df.rep:
            values = df[(df.phenotype == phenotype) & (df.transfer_n == t)]
            n = [list(values["rep"]).count(0), list(values["rep"]).count(
                1), list(values["rep"]).count(2), list(values["rep"]).count(3)]
            min_v.append(min(n))
            max_v.append(max(n))
            mean_v.append(np.mean(n))
        f_min.update({phenotype: np.array(min_v)/num_of_patients*100})
        f_max.update({phenotype: np.array(max_v)/num_of_patients*100})
        f_mean.update({phenotype: np.array(mean_v)/num_of_patients*100})
    return f_min, f_max, f_mean, T


def sort_legend_by(order, handles, labels):
    rank = list(range(7))
    order_dict = dict(zip(order, rank))
    tmp = pd.DataFrame()
    tmp["handles"] = handles
    tmp["labels"] = labels
    tmp["rank"] = tmp.apply(lambda x: order_dict[x["labels"]], axis=1)
    tmp.sort_values("rank", inplace=True)
    return list(tmp.handles), list(tmp.labels)


def plot_strategies(Data, strategies: list, color_dict, num_of_patients=94, sns_style=False, alpha=0.2, hatch=None, legend=False, handles=[[], []]):
    axs = handles[0]
    fig = handles[1]
    phenotypes = ["UI", "S", "A_r", "B_r", "A&B", "AB_r", "Fishy"]
    if len(axs) == 0:
        fig, axs = plt.subplots(2, 3, figsize=(
            25, 15), sharex=True, sharey=True, )
        fig.suptitle('TS-Plasmid: min/max plot', fontsize=30)
    j = 0
    k = 0
    for strategy in strategies:
        if sns_style:
            f_min, f_max, f_mean, T = get_min_max_for_phenotypes_sns(
                Data[Data.strategy == strategy], num_of_patients, phenotypes)
        else:
            f_min, f_max, f_mean, T = get_min_max_for_phenotypes(
                Data[Data.strategy == strategy], num_of_patients, phenotypes)
        if j == 2:
            j = 0
            k += 1
        for phenotype in Data.phenotype.unique():
            c = color_dict[phenotype]
            axs[j, k].fill_between(T, f_min[phenotype], f_max[phenotype],
                                   color=c, alpha=alpha, label=phenotype, hatch=hatch)
        for tick in axs[j, k].xaxis.get_major_ticks():
            tick.label.set_fontsize(18)
        for tick in axs[j, k].yaxis.get_major_ticks():
            tick.label.set_fontsize(18)
        if j == 1:
            axs[j, k].set_xlabel('transfer', fontsize=20)
        if k == 0:
            axs[j, k].set_ylabel(
                'fraction of patients with phenotype [%]', fontsize=20)
        axs[j, k].set_title(strategy, fontsize=25)
        j += 1

    handles, labels = axs[0, 0].get_legend_handles_labels()
    order = ["UI", "S", "A_r", "B_r", "A&B", "AB_r", "Fishy"]
    handles, labels = sort_legend_by(order, handles, labels)
    if legend:
        fig.legend(handles, labels, loc='lower center', ncol=len(
            labels), prop={'size': 20}, frameon=False)
    return axs, fig, handles, labels


def get_min_max_for_phenotypes_sns(df, num_of_patients, phenotypes):
    f_min = {}
    f_max = {}
    f_mean = {}
    T = df.transfer_n.unique()
    T.sort()
    for phenotype in phenotypes:
        min_v = []
        max_v = []
        mean_v = []
        for t in T:
            values = df[(df.phenotype == phenotype)
                        & (df.transfer_n == t)]["n"]
            min_v.append(min(values))
            max_v.append(max(values))
            mean_v.append(np.mean(values))
        f_min.update({phenotype: np.array(min_v)/num_of_patients*100})
        f_max.update({phenotype: np.array(max_v)/num_of_patients*100})
        f_mean.update({phenotype: np.array(mean_v)/num_of_patients*100})
    return f_min, f_max, f_mean, T


def plot_line_cloud(df, color_dict, alpha, num, ax, extreme=False):
    reps = list(df.rep.unique())
    selection = random.sample(reps, num)
    df = df[df.rep.isin(selection)]
    df.sort_values("transfer_n")
    for r in df.rep.unique():
        for pheno in df.phenotype.unique():
            line = df[(df.phenotype == pheno) & (df.rep == r)]
            if extreme:
                ax.plot(line["transfer_n"], line["f"], "-.",
                        linewidth=5, color=color_dict[pheno], alpha=alpha)
            else:
                ax.plot(line["transfer_n"], line["f"],
                        color=color_dict[pheno], alpha=alpha)


def plot_mean(ax, df_strat, color_dict, alpha, style, linewidth=1.5):
    cols = ["transfer_n", "phenotype", "f"]
    df_mean = df_strat[cols].groupby(
        cols[:-1]).mean().unstack(fill_value=0).stack().reset_index()
    for pheno in df_mean.phenotype.unique():
        df_pheno = df_mean[df_mean.phenotype == pheno]
        ax.plot(df_pheno["transfer_n"], df_pheno["f"], style,
                color=color_dict[pheno], alpha=alpha, linewidth=linewidth)
