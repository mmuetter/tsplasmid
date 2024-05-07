import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import string
import seaborn as sns


def make_cld(df, alpha=0.05):

    df["p-adj"] = df["p-adj"].astype(float)

    # Creating a list of the different treatment groups from Tukey's
    group1 = set(df.group1.tolist())  # Dropping duplicates by creating a set
    group2 = set(df.group2.tolist())  # Dropping duplicates by creating a set
    groupSet = group1 | group2  # Set operation that creates a union of 2 sets
    groups = sorted(list(groupSet))

    # Creating lists of letters that will be assigned to treatment groups
    letters = list(string.ascii_lowercase)[:len(groups)]
    cldgroups = letters

    # the following algoritm is a simplification of the classical cld,

    cld = pd.DataFrame(list(zip(groups, letters, cldgroups)))
    cld[3] = ""

    for row in df.itertuples():
        if df["p-adj"][row[0]] > (alpha):
            cld.iat[groups.index(df["group1"][row[0]]),
                    2] += cld.iat[groups.index(df["group2"][row[0]]), 1]
            cld.iat[groups.index(df["group2"][row[0]]),
                    2] += cld.iat[groups.index(df["group1"][row[0]]), 1]

        if df["p-adj"][row[0]] < (alpha):
            cld.iat[groups.index(df["group1"][row[0]]),
                    3] += cld.iat[groups.index(df["group2"][row[0]]), 1]
            cld.iat[groups.index(df["group2"][row[0]]),
                    3] += cld.iat[groups.index(df["group1"][row[0]]), 1]

    cld[2] = cld[2].apply(lambda x: "".join(sorted(x)))
    cld[3] = cld[3].apply(lambda x: "".join(sorted(x)))
    cld.rename(columns={0: "groups"}, inplace=True)

    # this part will reassign the final name to the group
    # for sure there are more elegant ways of doing this
    cld = cld.sort_values(cld.columns[2], key=lambda x: x.str.len())
    cld["labels"] = ""
    letters = list(string.ascii_lowercase)
    unique = []
    for item in cld[2]:

        for fitem in cld["labels"].unique():
            for c in range(0, len(fitem)):
                if not set(unique).issuperset(set(fitem[c])):
                    unique.append(fitem[c])
        g = len(unique)

        for kitem in cld[1]:
            if kitem in item:
                if cld["labels"].loc[cld[1] == kitem].iloc[0] == "":
                    cld["labels"].loc[cld[1] == kitem] += letters[g]

                # Checking if there are forbidden pairing (proposition of solution to the imperfect script)
                if kitem in ' '.join(cld[3][cld["labels"] == letters[g]]):
                    g = len(unique)+1

                # Checking if columns 1 & 2 of cld share at least 1 letter
                if len(set(cld["labels"].loc[cld[1] == kitem].iloc[0]).intersection(cld.loc[cld[2] == item, "labels"].iloc[0])) <= 0:
                    if letters[g] not in list(cld["labels"].loc[cld[1] == kitem].iloc[0]):
                        cld["labels"].loc[cld[1] == kitem] += letters[g]
                    if letters[g] not in list(cld["labels"].loc[cld[2] == item].iloc[0]):
                        cld["labels"].loc[cld[2] == item] += letters[g]

    cld = cld.sort_values("labels")
    # print(cld)
    # print('\n')
    cld.drop(columns=[1, 2, 3], inplace=True)
    # print(cld)
    # print('\n')
    # print('\n')
    return (cld.set_index("groups"))


def tukey(df, group_col, val_col, alpha=0.05):
    return pairwise_tukeyhsd(endog=df[val_col], groups=df[group_col], alpha=alpha)


def anova(df, group_col, var):
    groups = df[group_col].unique()
    V = []
    for group in groups:
        v = df.loc[df[group_col] == group, var].values
        V.append(v)
    stat = stats.f_oneway(*V)
    return stat


def summarize(Data, groups, label):
    summary = Data.groupby(groups)[label].agg(
        [(label+"_sum", 'sum'), (label+"_total", 'count'), (label+"_mean", 'mean')])
    return summary.reset_index()


def plot_tukey(summary, val_col, group_col, color, letters, name, lim=False, base_offset=0.04, ax=None, x_offset=0, alpha=0.05):
    means = summary.groupby(["strategy"]).mean()
    order = ['No treatment', 'Mono A', 'Mono B',
             'Combination', 'Cycling', 'Mixing']
    if bool(ax):
        sns.boxplot(showmeans=True,
                    meanline=True,
                    meanprops={'color': color, 'ls': '-', 'lw': 3},
                    medianprops={'visible': False},
                    whiskerprops={'visible': False},
                    x=group_col,
                    y=val_col,
                    data=summary,
                    showfliers=False,
                    showbox=False,
                    showcaps=False,
                    ax=ax,
                    order=order)
    else:
        # plot the mean line
        ax = sns.boxplot(showmeans=True,
                         meanline=True,
                         meanprops={'color': color, 'ls': '-', 'lw': 3},
                         medianprops={'visible': False},
                         whiskerprops={'visible': False},
                         x=group_col,
                         y=val_col,
                         data=summary,
                         showfliers=False,
                         showbox=False,
                         showcaps=False,
                         order=order)

    sns.stripplot(
        x=group_col,
        y=val_col,
        data=summary,
        color=color,
        s=15,
        alpha=alpha,
        ax=ax,
        order=order)

    if bool(lim):
        ax.set_ylim(*lim)
    else:
        ax.set_ylim(-0.1, 1.1)

    if type(letters) != type(None):
        for xtick in ax.get_xticks():
            label = ax.get_xticklabels()[xtick].get_text()
            pos_vec = summary.groupby(["strategy"]).mean()[val_col]-0.5
            if pos_vec[label]:
                offset = base_offset
            else:
                offset = -base_offset
            string = letters.loc[label, "labels"]
            ax.text(xtick+x_offset, means.loc[label, val_col]+offset, string,
                    horizontalalignment='center', color=color, weight='semibold', fontsize=15)

    return ax
