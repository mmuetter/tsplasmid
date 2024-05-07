from evaluation_scripts.base import rename_strings
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import levene
import numpy as np


class Model:
    def __init__(self, exp, models):
        exp.make_time_series()
        self.experiment = exp
        self.exp_data = exp.time_series.summary
        #  Models
        self.model_names = models.keys()
        self.model_types = list(models.values())
        self.models = []
        for model_name in models:
            self.load_model(model_name)

    def load_model(self, name):
        path = self.experiment.pathes["obj"]
        model = pd.read_pickle(os.path.join(path, name)).reset_index().rename(
            columns={"r": "rep", "x_hat_sim": "phenotype"})
        model = rename_strings(model, self.experiment.rename)
        self.models.append(model)

    def compare_models_to_exp(self):
        self.model_diff = []
        self.model_diff_arrays = []
        self.model_diff_mean = []
        l = len(self.models)
        self.model_diff_abs = []
        self.model_diff_var = []
        for i in range(0, l):
            self.exp_model_diff(i)
            diff_abs = self.model_diff[i]["diff_abs"]
            self.model_diff_mean.append(diff_abs.mean())
            self.model_diff_abs.append(diff_abs)
            self.model_diff_arrays.append(self.model_diff[i]["diff"])
            self.model_diff_var.append(np.var(self.model_diff[i]["diff"]))
        self.levene = levene(*self.model_diff_arrays)

    def exp_model_diff(self, i):
        model = self.models[i]
        exp_data = self.exp_data
        cols = ["transfer_n", "strategy", "phenotype", "f"]
        exp_grouped = exp_data[cols].groupby(cols[:-1]).mean()
        model_grouped = model[cols].groupby(cols[:-1]).mean()
        a = set(model_grouped.index)
        b = set(exp_grouped.index)

        differences = []
        for index in a.union(b):
            if (index in a) and (index in b):
                diff = exp_grouped.loc[index, "f"] - \
                    model_grouped.loc[index, "f"]
                row = dict(zip(["transfer_n", "strategy", "phenotype"], index))
                row.update({"diff": diff})
                differences.append(row)
        df = pd.DataFrame().from_records(differences)
        df["diff_abs"] = abs(df["diff"])
        self.model_diff.append(df)

    def plot_model_diff(self, figsize=(20, 10)):
        fig, axs = plt.subplots(1, len(self.models), figsize=figsize)
        for i, ax in enumerate(axs):
            sns.stripplot(data=self.model_diff[i], y="diff", ax=ax)
            ax.set_title(self.model_types[i])
