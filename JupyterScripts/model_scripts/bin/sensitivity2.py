import pandas as pd
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import matplotlib.patches as mpatches
import seaborn as sns
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os


class SensitivityAnalysis:
    def __init__(self, path, parameter_names=["turnover", "infection", "cA", "cB", "cS", "cU", "cAB"], name_add="_ab"):
        self.path = path
        self.file_path = os.path.join(
            path, "sensitivityAnalysis" + name_add + ".pkl")
        self.parameter_names = parameter_names

    def evaluate_parameter_sets(self, simulations, parameter):
        self.simulations = simulations
        self.strategies = simulations[list(simulations.keys())[
            0]].strategy.unique()

        self.parameter = parameter
        self.evaluations = []
        self.df = []
        for simulation, parameter in zip(list(self.simulations.values()), list(self.parameter.values())):
            parameter_eval = ParameterSetEvaluation(simulation, parameter)
            self.evaluations.append(parameter_eval)
            self.df.append(parameter_eval.summary_dict)
        self.data = pd.DataFrame().from_records(self.df)
        self.data = self.data.rename(columns={
            "U": "cU",
            "S": "cS",
            "A_r": "cA",
            "B_r": "cB",
            "AB_r": "cAB"
        })
        self.groups = self.data.winning_group.unique()
        colors = sns.color_palette('pastel', len(self.groups))
        self.group_colors = dict(zip(self.groups, colors))

    def lda(self):
        self.lda_data = self.data.copy()  # Create a copy of the data for LDA analysis
        X = self.lda_data[self.parameter_names]
        y = self.lda_data['winning_group']

        lda_model = LinearDiscriminantAnalysis(n_components=2)
        lda_data = lda_model.fit_transform(X, y)
        correlations = np.corrcoef(lda_data.T, X.T)

        self.lda_df = pd.DataFrame(
            lda_data, columns=['LDA Axis 1', 'LDA Axis 2'])
        self.lda_df['winning_group'] = y
        self.parameter_directions = pd.DataFrame(
            correlations[:2, 2:], columns=X.columns)

    def plot_lda(self):
        fig, ax = plt.subplots(figsize=(12, 12))
        for group in self.groups:
            strategy_data = self.lda_df[self.lda_df['winning_group'] == group]
            ax.scatter(strategy_data['LDA Axis 1'], strategy_data['LDA Axis 2'],
                       color=self.group_colors[group], label=group)
        ax.set_xlabel('LDA Axis 1')
        ax.set_ylabel('LDA Axis 2')
        ax.set_title('Preferable Strategy Map')
        ax.legend(title='Best Strategy')
        self.ax = ax
        self.fig = fig
        self.plot_original_axes()

    def plot_original_axes(self, arrow_length=2.5):
        for param in self.parameter_names:
            color = 'black'  # Set the color to black
            direction = self.parameter_directions[param].values
            self.ax.arrow(0, 0, direction[0] * arrow_length, direction[1] * arrow_length,
                          head_width=0.1, head_length=0.1, fc=color, ec=color)
            self.ax.text(direction[0] * arrow_length * 1.3, direction[1]
                         * arrow_length * 1.3, param, color=color, fontsize=15)

    def plot_lda_kde(self):
        fig, ax = plt.subplots(figsize=(12, 12))
        for group in self.groups:
            strategy_data = self.lda_df[self.lda_df['winning_group'] == group]
            x = strategy_data['LDA Axis 1']
            y = strategy_data['LDA Axis 2']

            # Plot the KDE using seaborn's kdeplot
            sns.kdeplot(x=x, y=y, ax=ax, label=group, fill=True)

        ax.set_xlabel('LDA Axis 1')
        ax.set_ylabel('LDA Axis 2')
        ax.set_title('Preferable Strategy Map (KDE)')

        # Create a custom legend
        legend_handles = []
        for group in self.groups:
            color = self.group_colors[group]
            patch = mpatches.Patch(color=color, label=group)
            legend_handles.append(patch)
        ax.legend(handles=legend_handles, title='Best Strategy')

        self.ax = ax
        self.fig = fig
        self.plot_original_axes()

    def save(self):
        with open(self.file_path, 'wb') as file:
            pickle.dump(self, file)

    def load(self):
        with open(self.file_path, 'rb') as file:
            obj = pickle.load(file)
        self.__dict__.update(obj.__dict__)


class ParameterSetEvaluation:
    def __init__(self, simulation, parameter_set):
        self.simulation = simulation
        self.parameter_set = parameter_set
        self.prepare_U_data()
        self.strategies = self.data["strategy"].unique()
        self.perform_anova()
        if self.anova_result.pvalue < 0.05:
            self.perform_tukey()
        self.label_winning_strategies()
        self.summerize()

    def prepare_U_data(self, exclude_noTreatment=True):
        self.t_end = t_end = max(self.simulation.transfer_n)
        sim_end = self.simulation[self.simulation.transfer_n.isin(
            [t_end, t_end-1, t_end-2, t_end-3])]
        sim_U = sim_end[sim_end.x_hat == "U"]
        if exclude_noTreatment:
            sim_U = sim_U[sim_U.strategy != "No treatment"]

        self.data = sim_U[["strategy", "f", "rep"]].groupby(
            ["strategy", "rep"]).mean().reset_index()

    def perform_anova(self):
        f_values = []
        for strategy in self.strategies:
            subset = self.data[self.data.strategy == strategy]
            f_value = subset["f"]
            f_values.append(f_value)
        self.anova_result = f_oneway(*f_values)

    def perform_tukey(self):
        self.tukey_results = tukey_result = pairwise_tukeyhsd(
            self.data["f"], self.data["strategy"])
        summary_df = pd.DataFrame(tukey_result.summary(
        ).data[1:], columns=tukey_result.summary().data[0])
        self.summary_df = summary_df

        if summary_df.reject.any():
            best_strategy = summary_df.iloc[summary_df['meandiff'].abs(
            ).idxmax()]['group1']
            self.best_strategy = best_strategy

            not_worse_strategies = summary_df[summary_df['group1']
                                              == best_strategy]
            significant_groups = not_worse_strategies[not_worse_strategies['reject'] == False]['group2'].tolist(
            )
            self.best_strategies = [best_strategy] + significant_groups

    def label_winning_strategies(self):
        if hasattr(self, 'best_strategies'):
            sorted_strategies = sorted(self.best_strategies)
            labeled_strategies = '-'.join(sorted_strategies)
            self.winning_group = labeled_strategies
        else:
            self.winning_group = "None"
            self.best_strategy = "None"

    def summerize(self):
        p = self.parameter_set.copy()
        self.summary_dict = p["community"]
        self.summary_dict.update({
            "turnover": p["turnover"],
            "infection": p["infection"],
            "winning_group": self.winning_group,
            "winning_strategy": self.best_strategy})
