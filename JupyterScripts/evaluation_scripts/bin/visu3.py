import pandas as pd
import seaborn as sns
import os
from itertools import product
from matplotlib import pyplot as plt


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


def load_model(name, pathes):
    model = pd.read_pickle(os.path.join(pathes["obj"], name)).reset_index().rename(
        columns={"r": "rep", "x_hat_sim": "phenotype"})
    return model


def plot_exp_data(ax, df, color_dict, scale=2):
    r = df.rep.unique()
    p = df.phenotype.unique()
    s = df.strategy.unique()[0]
    df = df.dropna()  # drop rows with NaN values

    L = []
    for t in range(0, max(df.transfer_n+1)):
        if t not in df.transfer_n.unique():
            for phenotype, rep in product(p, r):
                L.append({"rep": rep, "phenotype": phenotype,
                         "f": None, "strategy": s, "transfer_n": t})
    fill = pd.DataFrame(L)
    df = pd.concat([fill, df], sort=True)
    print("1")
    sns.pointplot(data=df,
                  ax=ax,
                  x="transfer_n",
                  y="f",
                  scale=scale,
                  hue="phenotype",
                  errorbar=("pi", 100), join=False, dodge=0.2,
                  palette=color_dict,
                  )


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
