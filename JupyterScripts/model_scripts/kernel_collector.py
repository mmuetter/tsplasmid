import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import pickle
import os


class GridPoint:
    def __init__(self, x, y, i, j, r, data, strategy_columns):
        self.x = x
        self.y = y
        self.i = i
        self.j = j
        self.r = r
        self.data = data
        self.strategy_columns = strategy_columns[:5]
        self.collect_values()
        if not self.empty:
            self.anova()
            if self.significance:
                self.tukey()
            else:
                self.winning_group = "No winner"
                self.winners = []

    def collect_values(self):
        lda_axis_1 = self.data["LDA_axis_1"]
        lda_axis_2 = self.data["LDA_axis_2"]
        f_values = self.data.drop(columns=["LDA_axis_1", "LDA_axis_2"])
        distances = distance.cdist(
            [[self.x, self.y]], np.column_stack((lda_axis_1, lda_axis_2)))
        indices = np.where(distances <= self.r)[1]
        self.collected_values = f_values.iloc[indices]
        self.df = pd.DataFrame(self.collected_values,
                               columns=self.collected_values.columns)
        self.num_of_points = len(indices)
        self.empty = self.df.empty
        if self.empty:
            self.num_of_points = 0

    def anova(self):
        anova_data = self.df[self.strategy_columns[:5]]
        _, p_value = f_oneway(*anova_data.values.T)
        self.significance = p_value < 0.05

    def tukey(self):
        strategy_data = [self.df[col] for col in self.strategy_columns[:5]]
        tukey_results = pairwise_tukeyhsd(np.concatenate(
            strategy_data), np.repeat(self.strategy_columns[:5], len(self.df)))
        best_strategy = tukey_results.groups[np.argmin(
            tukey_results.meandiffs)]
        self.winners = [best_strategy]
        for row in tukey_results.summary().data:
            group1, group2, reject = row[0], row[1], row[-1]
            if group1 == best_strategy and not reject:
                self.winners.append(group2)

        self.winners.sort()  # Sort the winners alphabetically
        self.winning_group = "_".join(self.winners)


class KernelCollector:
    def __init__(self, sensitivity, path, r=0.1, name_add=""):
        self.sensitivity = sensitivity
        self.r = r
        self.strategy_columns = sensitivity.strategies
        self.parameter_directions = sensitivity.parameter_directions
        self.path = path
        self.file_path = os.path.join(
            path, "kernels" + name_add + ".pkl")

    def create_grid_points(self, num_ticks):
        self.num_ticks = num_ticks
        lda_axis_1 = self.sensitivity.evaluations_df["LDA_axis_1"]
        lda_axis_2 = self.sensitivity.evaluations_df["LDA_axis_2"]
        min_x, max_x = min(lda_axis_1), max(lda_axis_1)
        min_y, max_y = min(lda_axis_2), max(lda_axis_2)
        self.x_vals = np.linspace(min_x, max_x, self.num_ticks)
        self.y_vals = np.linspace(min_y, max_y, self.num_ticks)
        self.grid_points = []
        for i, x in enumerate(self.x_vals):
            print(i)
            for j, y in enumerate(self.y_vals):
                grid_point = GridPoint(
                    x, y, i, j, self.r, self.sensitivity.evaluations_df, self.strategy_columns)
                self.grid_points.append(grid_point)

    def save(self):
        with open(self.file_path, 'wb') as file:
            pickle.dump(self, file)

    def load(self):
        with open(self.file_path, 'rb') as file:
            obj = pickle.load(file)
        self.__dict__.update(obj.__dict__)
